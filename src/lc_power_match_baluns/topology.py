# This file is part of lc-power-match-baluns.
# Copyright © 2023 Technical University of Denmark (developed by Rasmus Jepsen)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

"""This module encapsulates the functionality of LC-balun topologies.
"""

import lcapy
import sympy
import lc_power_match_baluns.oneport
from collections import abc
import math
import numpy as np
import skrf
from multimethod import multimethod

class BalunTopology:
  """A power matching LC-balun topology"""

  name: str = ""

  netlist: str = ""

  num_elements: int = 0

  _two_port_z_params_func = None
  _three_port_z_params_func = None

  @multimethod
  def calculate_two_port_impedance_parameters(cls, # pylint: disable=no-self-argument
      element_impedances: abc.Sequence[complex]):
    """Calculates two-port impedance parameters from element impedances

    Args:
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms

    Returns:
        The z parameters in ohms
    """
    element_impedances_array = np.array(element_impedances, dtype=complex)
    if cls._two_port_z_params_func is None:
      num_elements = len(element_impedances)
      circuit = cls.lcapy_circuit()
      impedance_symbols = [lcapy.symbol(f"Z{i + 1}").sympy for i in range(num_elements)]
      z_params = circuit.Zparamsn("1_0", "3_0", "2_0", 0).sympy
      cls._two_port_z_params_func = sympy.lambdify([tuple(impedance_symbols)], z_params)
    result = cls._two_port_z_params_func(element_impedances_array) # pylint: disable=not-callable
    return result # type: ignore [return-value]
  
  @classmethod
  @multimethod
  def calculate_two_port_impedance_parameters(cls, # pylint: disable=function-redefined
      element_impedances: abc.Sequence[abc.Sequence[complex]]):
    """Calculates two-port impedance parameters from element impedances

    Args:
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency

    Returns:
        The z parameters in ohms over frequency
    """
    return np.array([cls.calculate_two_port_impedance_parameters(impedances) for impedances in zip(*element_impedances)]) # pylint: disable=no-value-for-parameter

  @multimethod
  def calculate_two_port_scattering_parameters(cls, # pylint: disable=no-self-argument
      zb: complex, zu: complex, element_impedances: abc.Sequence[complex],
      s_def="power"):
    """Calculates two-port scattering parameters from element impedances

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Returns:
        The s parameters
    """
    s_params = cls.calculate_two_port_scattering_parameters([zb], [zu], list(zip(element_impedances)), s_def)
    return s_params[0, :, :]
  
  @classmethod
  @multimethod
  def calculate_two_port_scattering_parameters(cls, # pylint: disable=function-redefined
      zb: abc.Sequence[complex], zu: abc.Sequence[complex], element_impedances: abc.Sequence[abc.Sequence[complex]],
      s_def="power"):
    """Calculates two-port scattering parameters from element impedances

    Args:
        zb (abc.Sequence[complex]): Impedance presented to the balanced port in ohms over frequency
        zu (abc.Sequence[complex]): Impedance presented to the unbalanced port in ohms over frequency
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Note: The frequencies for zb, zu and element_impedances are assumed to be the same

    Returns:
        The s parameters over frequency
    """
    twoport_z = cls.calculate_two_port_impedance_parameters(element_impedances) # pylint: disable=no-value-for-parameter
    network = skrf.Network.from_z(twoport_z, s_def = s_def, z0=list(zip(zb, zu)))
    result = network.s
    return result # type: ignore [return-value]

  @multimethod
  def calculate_insertion_loss(cls, # pylint: disable=no-self-argument
      zb: complex, zu: complex, element_impedances: abc.Sequence[complex],
      s_def="power") -> float:
    """Calculates differential mode insertion loss from element impedances

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Returns:
        The insertion loss in decibels
    """
    insertion_loss = cls.calculate_insertion_loss([zb], [zu], list(zip(element_impedances)), s_def)
    return insertion_loss[0]
  
  @classmethod
  @multimethod
  def calculate_insertion_loss(cls, # pylint: disable=function-redefined
      zb: abc.Sequence[complex], zu: abc.Sequence[complex],
      element_impedances: abc.Sequence[abc.Sequence[complex]], s_def="power") -> abc.Sequence[float]:
    """Calculates differential mode insertion loss from element impedances

    Args:
        zb (abc.Sequence[complex]): Impedance presented to the balanced port in ohms over frequency
        zu (abc.Sequence[complex]): Impedance presented to the unbalanced port in ohms over frequency
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".
    
    Note: The frequencies for zb, zu and element_impedances are assumed to be the same

    Returns:
        The insertion loss in decibels over frequency
    """
    twoport_s = cls.calculate_two_port_scattering_parameters(zb, zu, element_impedances, s_def)
    s21 = twoport_s[:, 1, 0]
    insertion_loss = -20 * np.log10(np.abs(s21))
    return insertion_loss

  @multimethod
  def calculate_return_losses(cls, # pylint: disable=no-self-argument
      zb: complex, zu: complex, element_impedances: abc.Sequence[complex],
      s_def="power") -> tuple[float, float]:
    """Calculates differential mode return losses from element impedances

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Returns:
        The return losses in decibels (balanced port, unbalanced port)
    """
    balanced_return_loss, unbalanced_return_loss = cls.calculate_return_losses([zb], [zu], list(zip(element_impedances)), s_def)
    return (balanced_return_loss[0], unbalanced_return_loss[0])
  
  @classmethod
  @multimethod
  def calculate_return_losses(cls, # pylint: disable=function-redefined
      zb: abc.Sequence[complex], zu: abc.Sequence[complex],
      element_impedances: abc.Sequence[abc.Sequence[complex]], s_def="power") -> tuple[abc.Sequence[float], abc.Sequence[float]]:
    """Calculates differential mode return losses from element impedances

    Args:
        zb (abc.Sequence[complex]): Impedance presented to the balanced port in ohms over frequency
        zu (abc.Sequence[complex]): Impedance presented to the unbalanced port in ohms over frequency
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Note: The frequencies for zb, zu and element_impedances are assumed to be the same

    Returns:
        The return losses in decibels (balanced port, unbalanced port) over frequency
    """
    twoport_s = cls.calculate_two_port_scattering_parameters(zb, zu, element_impedances, s_def)
    s11 = twoport_s[:, 0, 0]
    s22 = twoport_s[:, 1, 1]
    balanced_return_loss = -20 * np.log10(np.abs(s11))
    unbalanced_return_loss = -20 * np.log10(np.abs(s22))
    return (balanced_return_loss, unbalanced_return_loss)

  @multimethod
  def calculate_three_port_impedance_parameters(cls, # pylint: disable=no-self-argument
      element_impedances: abc.Sequence[complex]):
    """Calculates three-port impedance parameters from element impedances

    Args:
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms

    Raises:
        NotImplementedError: If the topology has not implemented this method

    Returns:
        The z parameters in ohms
    """
    if cls._three_port_z_params_func is None:
      num_elements = len(element_impedances)
      circuit = cls.lcapy_circuit()
      impedance_symbols = [lcapy.symbol(f"Z{i + 1}").sympy for i in range(num_elements)]
      z_params = circuit.Zparamsn("2_0", 0, "1_0", 0, "3_0", 0).sympy
      cls._three_port_z_params_func = sympy.lambdify([tuple(impedance_symbols)], z_params)
    result = cls._three_port_z_params_func(element_impedances) # pylint: disable=not-callable
    return result # type: ignore [return-value]
  
  @classmethod
  @multimethod
  def calculate_three_port_impedance_parameters(cls, # pylint: disable=function-redefined
      element_impedances: abc.Sequence[abc.Sequence[complex]]):
    """Calculates three-port impedance parameters from element impedances

    Args:
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency

    Raises:
        NotImplementedError: If the topology has not implemented this method

    Returns:
        The z parameters in ohms over frequency
    """
    return np.array([cls.calculate_three_port_impedance_parameters(impedances) for impedances in zip(*element_impedances)]) # pylint: disable=no-value-for-parameter
  
  @multimethod
  def calculate_three_port_scattering_parameters(cls, # pylint: disable=no-self-argument
      zb: complex, zu: complex, element_impedances: abc.Sequence[complex],
      s_def="power"):
    """Calculates three-port scattering parameters from element impedances

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms
        element_impedances (abc.Sequence[complex]): The impedances of each element
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Returns:
        The s parameters
    """
    s_params = cls.calculate_three_port_scattering_parameters([zb], [zu], list(zip(element_impedances)), s_def)
    return s_params[0, :, :]
  
  @classmethod
  @multimethod
  def calculate_three_port_scattering_parameters(cls, # pylint: disable=function-redefined
      zb: abc.Sequence[complex], zu: abc.Sequence[complex],
      element_impedances: abc.Sequence[abc.Sequence[complex]], s_def="power"):
    """Calculates three-port scattering parameters from element impedances

    Args:
        zb (abc.Sequence[complex]): Impedance presented to the balanced port in ohms over frequency
        zu (abc.Sequence[complex]): Impedance presented to the unbalanced port in ohms over frequency
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Note: The frequencies for zb, zu and element_impedances are assumed to be the same

    Returns:
        The s parameters over frequency
    """
    threeport_z = cls.calculate_three_port_impedance_parameters(element_impedances) # pylint: disable=no-value-for-parameter
    network = skrf.Network.from_z(threeport_z, s_def = s_def, z0=list(zip(zu, np.array(zb) / 2, np.array(zb) / 2)))
    result = network.s
    return result # type: ignore [return-value]
  
  @multimethod
  def calculate_cmrr(cls, # pylint: disable=no-self-argument
      zb: complex, zu: complex, element_impedances: abc.Sequence[complex],
      s_def="power") -> float:
    """Calculates common mode rejection ratio from element impedances

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms
        element_impedances (abc.Sequence[complex]): The impedances of each element in ohms
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Returns:
        The common mode rejection ratio in decibels
    """
    cmrr_db = cls.calculate_cmrr([zb], [zu], list(zip(element_impedances)), s_def)
    return cmrr_db[0]
  
  @classmethod
  @multimethod
  def calculate_cmrr(cls, # pylint: disable=function-redefined
      zb: abc.Sequence[complex], zu: abc.Sequence[complex],
      element_impedances: abc.Sequence[abc.Sequence[complex]], s_def="power") -> abc.Sequence[float]:
    """Calculates common mode rejection ratio from element impedances

    Args:
        zb (abc.Sequence[complex]): Impedance presented to the balanced port in ohms over frequency
        zu (abc.Sequence[complex]): Impedance presented to the unbalanced port in ohms over frequency
        element_impedances (abc.Sequence[abc.Sequence[complex]]): A sequence of all elements,
            with the inner sequence representing the impedance of an element in ohms over frequency
        s_def (str, optional): The scattering parameter definition to use from scikit-rf. Defaults to "power".

    Note: The frequencies for zb, zu and element_impedances are assumed to be the same

    Returns:
        The common mode rejection ratio in decibels over frequency
    """
    threeport_s = cls.calculate_three_port_scattering_parameters(zb, zu, element_impedances, s_def)
    linear_cmrr = (threeport_s[:, 1, 0] - threeport_s[:, 2, 0]) / (threeport_s[:, 1, 0] + threeport_s[:, 2, 0])
    cmrr_db = 20 * np.log10(np.abs(linear_cmrr))
    return cmrr_db

  @classmethod
  def lcapy_circuit(cls) -> lcapy.Circuit:
    """Creates an Lcapy circuit form of this topology from its netlist.

    The circuit form uses generalized impedance elements for each element.

    Returns:
        lcapy.Circuit: The Lcapy circuit form
    """
    return lcapy.Circuit(cls.netlist)

  @classmethod
  def lc_netlist(cls, components: abc.Sequence[lc_power_match_baluns.oneport.SimpleLosslessOnePort]) -> str:
    """Creates a netlist of this topology with inductors and capacitors.

    The network uses an inductor or a capacitor for each element depending on the values of each component.

    Args:
        components (abc.Sequence[five_element.oneport.SimpleLosslessOnePort]): Components used to determine how each element is implemented.

    Returns:
        str: The netlist with LC elements
    """
    netlist = cls.netlist
    for component in components:
      netlist = netlist.replace(f"Z{component.index}", component.symbol)
    return netlist
  
  @classmethod
  def lcapy_lc_circuit(cls, components: abc.Sequence[lc_power_match_baluns.oneport.SimpleLosslessOnePort]) -> lcapy.Circuit:
    """Creates an Lcapy circuit form of this topology from its netlist.

    The network uses an inductor or a capacitor for each element depending on the values of each component.

    Args:
        components (abc.Sequence[five_element.oneport.SimpleLosslessOnePort]): Components used to determine how each element is implemented.

    Returns:
        lcapy.Circuit: The Lcapy circuit form with LC elements
    """
    return lcapy.Circuit(cls.lc_netlist(components))
  
  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    raise NotImplementedError

  @classmethod
  def calculate_elements_from_impedances(cls, zb: complex, zu: complex) -> abc.Sequence[tuple[float, ...]]:
    """Calculates element reactances from impedances at each port for power matching and common-mode rejection

    Args:
        zb (complex): Impedance presented to the balanced port in ohms
        zu (complex): Impedance presented to the unbalanced port in ohms

    Raises:
        NotImplementedError: If the topology has not implemented this method

    Returns:
        Sequence[tuple[float, ...]]: A sequence containing the element reactances in ohms for each solution
    """
    rb, xb, ru, xu = zb.real, zb.imag, zu.real, zu.imag
    return cls._calculate_elements_from_impedances(rb, xb, ru, xu)

class ExtendedTTopology(BalunTopology):
  """Extended T balun topology for matching complex impedances.
  Fritz et al. [1] have previously considered an equivalent topology and its design equations.
  The design equations that have been used are those derived in the derivations folder independently.

  [1] M. Fritz, M. Handtmann, and P. Bradley, "Four lc element balun," English, pat. 9 106 204, 2013. [Online]. Available: https://patents.google.com/patent/US9106204.
  """

  name = "Extended T"

  netlist = """
            Z1 1_0 5_0; right
            Z2 5_0 3_1; down
            W 3_0 3_1; right
            Z4 3_1 0; right
            Z3 5_0 2_0; right
            W 0 0_0; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 4

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zb = math.sqrt(rb ** 2 + xb ** 2)
    x1_1 = -zb * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    x1_2 = zb * math.sqrt(ru / rb)
    x2_1 = zb * math.sqrt(ru / rb)
    x2_2 = -zb * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    x3_1 = ru * xb / rb - xu - zb / 2 * math.sqrt(ru / rb)
    x3_2 = ru * xb / rb - xu + zb / 2 * math.sqrt(ru / rb)
    x4_1 = -zb / 2 * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    x4_2 = zb / 2 * math.sqrt(ru / rb)
    return [(x1_1, x2_1, x3_1, x4_1), (x1_2, x2_2, x3_2, x4_2)]
  
class ExtendedPiTopology(BalunTopology):
  """Extended Pi balun topology for matching complex impedances.
  Bradley and Frank [2] have previously considered a special case of this topology.
  The design equations that have been used are those derived in the derivations folder independently.

  [2] P. Bradley and M. Frank, "Combined balun and impedance matching circuit," English, pat. 8 633 781, 2010. [Online]. Available: https://patents.google.com/patent/US8633781B2/en.
  """

  name = "Extended Pi"

  netlist = """
            Z1 1_0 3_0; down
            Z2 1_0 2_0; right
            Z3 2_0 3_1; down
            Z4 3_1 0; right
            W 3_0 3_1; right
            W 2_0 2_1; right
            W 1_1 1_0; right
            W 3_2 3_0; right
            W 0 0_0; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 4

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zb = math.sqrt(rb ** 2 + xb ** 2)
    x1_1 = float(np.float64(2 * zb ** 2 * ru) / (2 * xu * rb - 2 * ru * xb - zb * math.sqrt(ru * rb)))
    x1_2 = float(np.float64(2 * zb ** 2 * ru) / (2 * xu * rb - 2 * ru * xb + zb * math.sqrt(ru * rb)))
    x2_1 = zb * math.sqrt(ru / rb)
    x2_2 = -zb * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    x3_1 = -zb * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    x3_2 = zb * math.sqrt(ru / rb)
    x4_1 = zb / 2 * math.sqrt(ru / rb)
    x4_2 = -zb / 2 * math.sqrt(ru / rb) # pylint: disable=invalid-unary-operand-type
    return [(x1_1, x2_1, x3_1, x4_1), (x1_2, x2_2, x3_2, x4_2)]

class LatticeTopology(BalunTopology):
  """Generalized lattice balun topology for matching complex impedances.
  Symmetric lattice baluns have previously been used for real-real matching [3],
  though they require extra elements if used for complex impedance matching.
  This topology relaxes the symmetry constraint to provide complex impedance matching with four elements.
  Apel and Page [4] have also previously considered a balun network that resembles a specific case of an asymmetric lattice balun.

  The design equations that have been used are those derived in the derivations folder independently.

  [3] C Lorenz AG, "Circuit arrangement for the transition from a symmetrical electrical arrangement to an asymmetrical one, in particular in the case of high-frequency arrangements," Germany Patent 603 816, April 1, 1932. [Online]. Available: https://patents.google.com/patent/DE603816C/en
  [4] T. R. Apel and C. E. Page, "Lumped parameter balun," English, pat. 5 574 411, 1995. [Online]. Available: https://patents.google.com/patent/US5574411A/en.
  """

  name = "Lattice"

  netlist = """
            Z2 5_1 6_0; right
            W 6_0 6_1; right
            Z4 6_1 4_2; rotate=225
            W 4_2 4_1; rotate=225
            W 4_1 4_0; right
            Z3 4_0 0_3; right
            Z1 5_2 0_0; rotate=-45
            W 5_2 5_1; right
            W 0_1 0_3; rotate=-45
            W 0_0 0_1; right
            W 0_3 0_2; right
            W 0_2 0; right
            W 1_0 5_2; right
            W 3_0 4_1; right
            W 6_1 2_0; right
            W 0 0_4; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 4

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zb = math.sqrt(rb ** 2 + xb ** 2)
    x1 = float(np.float64(ru * zb ** 2) / (2 * rb * xu - 2 * ru * xb - zb * math.sqrt(ru * rb)))
    x2 = zb * math.sqrt(ru / rb)
    x3 = float(np.float64(ru * zb ** 2) / (2 * rb * xu - 2 * ru * xb + zb * math.sqrt(rb * ru)))
    x4 = -zb * math.sqrt(ru / rb)
    return [(x1, x2, x3, x4)]

class DipperTopology(BalunTopology):
  name = "Dipper"

  netlist = """
Z2 5 4; down
W 5 2; right
Z4 2 0_2; down
Z3 4 0_2; right
W 3_0 4; right
Z1 1_0 5; right
W 2 2_0; right
W 0_2 0; right
W 0 0_1; down=0.1, ground
; label_nodes=none
"""

  num_elements = 4
  
  @classmethod
  def _calc_dipper_reactances_special(cls, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zu2 = ru ** 2 + xu ** 2
    x1 = -xb / 2 + 2 * xu
    x2 = xb / 2 - 2 * xu
    x3 = -xb / 4 + xu
    if x1 == 0 or x2 == 0 or x3 == 0:
      return []
    x4 = float(np.float64(zu2) * (4 * xu - xb) / (4 * ru ** 2 - 4 * xu ** 2 + 2 * xb * xu))
    return [(x1, x2, x3, x4)]

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zu2 = ru ** 2 + xu ** 2
    zb2 = rb ** 2 + xb ** 2
    delta = 4 * zu2 - rb * ru
    if delta < 0:
      return []
    factor = rb - 4 * ru
    if factor == 0:
      return cls._calc_dipper_reactances_special(xb, ru, xu)
    x1_1 = -xb / 2 - math.sqrt(rb * delta / ru) / 2
    x1_2 = -xb / 2 + math.sqrt(rb * delta / ru) / 2
    x2_1 = xb / 2 + math.sqrt(rb * delta / ru) / 2
    x2_2 = xb / 2 - math.sqrt(rb * delta / ru) / 2
    x3_1 = -xb / 4 - math.sqrt(rb * delta / ru) / 4
    x3_2 = -xb / 4 + math.sqrt(rb * delta / ru) / 4
    x4_1 = float(np.float64(ru * zb2 / rb ** 2 * factor) / (xb + 4 * xu - 4 * ru * xb / rb + math.sqrt(rb * delta / ru)) - xu - ru * xb / rb)
    x4_2 = float(np.float64(ru * zb2 / rb ** 2 * factor) / (xb + 4 * xu - 4 * ru * xb / rb - math.sqrt(rb * delta / ru)) - xu - ru * xb / rb)
    return [(x1_1, x2_1, x3_1, x4_1), (x1_2, x2_2, x3_2, x4_2)]

class YuTopology(BalunTopology):
  name = "Yu"

  netlist = """
            Z1 1_0 2_1; right
            Z3 2_1 4_0; down
            Z2 3_0 4_0; right
            Z4 4_0 0; right
            W 2_1 2_0; right
            W 0 0_0; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 4

  @classmethod
  def _calc_yu_reactances_special(cls, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    if xu == 0:
      return []
    zu2 = ru ** 2 + xu ** 2
    x1 = 2 * xu - xb / 2
    x2 = xu - ru ** 2 / xu - xb / 2
    x3 = -xu - ru ** 2 / xu
    x4 = zu2 / 2 / xu
    return [(x1, x2, x3, x4)]
  
  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zu2 = ru ** 2 + xu ** 2
    delta = 4 * zu2 - rb * ru
    if delta < 0:
      return []
    denominator = 4 * ru - rb
    if denominator == 0:
      return cls._calc_yu_reactances_special(xb, ru, xu)
    x1_1 = -xb / 2 + math.sqrt(rb * delta / ru) / 2
    x1_2 = -xb / 2 - math.sqrt(rb * delta / ru) / 2
    x2_1 = 2 * rb * xu / (4 * ru - rb) - xb / 2 - math.sqrt(rb ** 3 * delta / ru) / 2 / denominator
    x2_2 = 2 * rb * xu / (4 * ru - rb) - xb / 2 + math.sqrt(rb ** 3 * delta / ru) / 2 / denominator
    x3_1 = (2 * rb * xu - 2 * math.sqrt(ru * rb * delta)) / denominator
    x3_2 = (2 * rb * xu + 2 * math.sqrt(ru * rb * delta)) / denominator
    x4_1 = -(rb * xu - math.sqrt(ru * rb * delta)) / denominator
    x4_2 = -(rb * xu + math.sqrt(ru * rb * delta)) / denominator
    return [(x1_1, x2_1, x3_1, x4_1), (x1_2, x2_2, x3_2, x4_2)]
  
class ReverseYuTopology(BalunTopology):
  name = "Reverse Yu"

  netlist = """
            W 1_0 1_1; right
            Z2 1_1 4_0; down
            Z1 3_0 4_0; right
            Z4 4_0 0; right
            Z3 1_1 2_0; right
            W 0 0_0; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 4

  @classmethod
  def _calc_reverse_yu_reactances_special(cls, xb: float, ru: float, xu: float):
    if xb == 0:
      return []
    x1 = -4 * ru ** 2 / xb - xb / 4
    x2 = -4 * ru ** 2 / xb - xb / 4
    x3 = xb / 4 - xu
    x4 = 2 * ru ** 2 / xb + xb / 8
    return [(x1, x2, x3, x4)]

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    zb2 = rb ** 2 + xb ** 2
    delta = -4 * ru * rb + zb2
    if delta < 0:
      return []
    denominator = 4 * ru - rb
    if denominator == 0:
      return cls._calc_reverse_yu_reactances_special(xb, ru, xu)
    x1_1 = (-2 * ru * xb - math.sqrt(ru * rb * delta)) / denominator
    x1_2 = (-2 * ru * xb + math.sqrt(ru * rb * delta)) / denominator
    x2_1 = (-2 * ru * xb - math.sqrt(ru * rb * delta)) / denominator
    x2_2 = (-2 * ru * xb + math.sqrt(ru * rb * delta)) / denominator
    x3_1 = -xu - math.sqrt(ru * delta / rb) / 2
    x3_2 = -xu + math.sqrt(ru * delta / rb) / 2
    x4_1 = (2 * ru * xb + math.sqrt(ru * rb * delta)) / 2 / denominator
    x4_2 = (2 * ru * xb - math.sqrt(ru * rb * delta)) / 2 / denominator
    return [(x1_1, x2_1, x3_1, x4_1), (x1_2, x2_2, x3_2, x4_2)]
  
class TraditionalLatticeTopology(BalunTopology):
  """The lattice topology from [3] extended with three elements to enable it to power match complex impedances.

  [3] C Lorenz AG, "Circuit arrangement for the transition from a symmetrical electrical arrangement to an asymmetrical one, in particular in the case of high-frequency arrangements," Germany Patent 603 816, April 1, 1932. [Online]. Available: https://patents.google.com/patent/DE603816C/en
  """

  name = "Traditional Lattice"

  netlist = """
            Z1 1_0 5_0; right
            Z4 5_1 6_0; right
            W 6_0 6_1; right
            Z6 6_1 4_2; rotate=225
            W 4_2 4_1; rotate=225
            W 4_1 4_0; right
            Z5 4_0 0_3; right
            Z3 5_2 0_0; rotate=-45
            W 5_2 5_1; right
            W 0_1 0_3; rotate=-45
            W 0_0 0_1; right
            W 0_3 0_2; right
            W 0_2 0; right
            W 5_0 5_2; right
            Z2 3_0 4_1; right
            Z7 6_1 2_0; right
            W 0 0_4; down=0.1, ground
            ; label_nodes=none
        """

  num_elements = 7

  @classmethod
  def _calculate_elements_from_impedances(cls, rb: float, xb: float, ru: float, xu: float) -> abc.Sequence[tuple[float, ...]]:
    x1 = x2 = -xb / 2
    x4 = x5 = math.sqrt(rb * ru)
    x3 = x6 = -math.sqrt(rb * ru) # pylint: disable=invalid-unary-operand-type
    x7 = -xu
    return [(x1, x2, x3, x4, x5, x6, x7)]