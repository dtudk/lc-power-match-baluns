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

"""This module encapsulates the behaviour of simple lossless one-port elements.
"""

from __future__ import annotations
import math
from collections import abc

class SimpleLosslessOnePort:
  """Represents either an inductor or a capacitor"""

  prefix: str = ""

  index: int | None = None

  value: float = 0

  unit: str = ""

  @property
  def symbol(self) -> str:
    """Symbol for the element 

    Returns:
        str: _description_
    """
    return f"{self.prefix}{self.index}"
  
  def __init__(self, value: float, index: int | None = None):
    """Creates a simple lossless one-port with a given value and index.

    Args:
        value (float): The component value for this element.
        index (int | None, optional): The index for this element. Defaults to None.
    """
    super().__init__()
    self.index = index
    self.value = value
  
  def __repr__(self) -> str:
    return f"{self.symbol} = {self.value} {self.unit}"
  
  def calc_reactance(self, frequency: float) -> float:
    """Calculates the reactance of this element at a given frequency

    Args:
        frequency (float): The frequency to calculate the reactance at in Hertz

    Raises:
        NotImplementedError: If the element has not implemented this method

    Returns:
        float: The element reactance in Ohms
    """
    raise NotImplementedError

  @classmethod
  def from_reactance_at_frequency(cls, reactance: float, frequency: float, index: int) -> SimpleLosslessOnePort:
    """Creates an element instance with a given reactance at a given frequency

    Creates a capacitor if the reactance is negative, creates an inductor otherwise

    Args:
        reactance (float): The reactance of the element in Ohms
        frequency (float): The frequency in Hertz
        index (int): The index of the created element

    Returns:
        SimpleLosslessOnePort: The created element
    """
    if reactance < 0:
      return Capacitor.from_reactance_at_frequency(reactance, frequency, index)
    else:
      return Inductor.from_reactance_at_frequency(reactance, frequency, index)
  
  @classmethod
  def from_reactances_at_frequency(cls, reactances: abc.Sequence[float], frequency: float) -> list[SimpleLosslessOnePort]:
    """Creates a sequence of elements from a sequence of reactances at a given frequency

    Each element is created as an inductor or capacitor depending on the sign of each reactance

    Args:
        reactances (abc.Sequence[float]): The sequence of reactances in Ohms
        frequency (float): The frequency in Hertz

    Returns:
        list[SimpleLosslessOnePort]: A list of the created elements
    """
    return [SimpleLosslessOnePort.from_reactance_at_frequency(reactance, frequency, i + 1) for i, reactance in enumerate(reactances)]

class Capacitor(SimpleLosslessOnePort):
  prefix: str = "C"

  unit: str = "F"

  def calc_reactance(self, frequency: float) -> float:
    return -1 / (2 * math.pi * frequency * self.value)
  
  @classmethod
  def from_reactance_at_frequency(cls, reactance: float,
      frequency: float, index: int | None = None) -> Capacitor:
    value = -1 / (2 * math.pi * frequency * reactance)
    return Capacitor(value, index)

class Inductor(SimpleLosslessOnePort):
  prefix: str = "L"

  unit: str = "H"
  
  def calc_reactance(self, frequency: float) -> float:
    return 2 * math.pi * frequency * self.value
  
  @classmethod
  def from_reactance_at_frequency(cls, reactance: float,
      frequency: float, index: int | None = None) -> Inductor:
    value = reactance / (2 * math.pi * frequency)
    return Inductor(value, index)