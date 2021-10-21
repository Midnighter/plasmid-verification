# Copyright (c) 2021, Moritz E. Beber.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""Provide a container class for a local sequence alignment."""


import enum
import re
from dataclasses import dataclass
from typing import ClassVar, Iterable, Pattern, Tuple

from Bio.Seq import Seq


class CIGAROperation(enum.Enum):
    """
    Define CIGAR operations.

    See https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format for more info.

    """

    match = "M"
    insertion = "I"
    deletion = "D"
    skipped = "N"


@dataclass(frozen=True)
class SequenceAlignment:

    score: float
    query_sequence: Seq
    query_begin: int
    query_end: int
    target_sequence: Seq
    target_begin: int
    target_end: int
    cigar: str
    aligned_query_sequence: Seq
    aligned_target_sequence: Seq

    _cigar_pattern: ClassVar[Pattern] = re.compile(r"(?P<num>\d+)(?P<code>[MNDI])")

    def iter_cigar(self) -> Iterable[Tuple[int, CIGAROperation]]:
        """"""
        return (
            (int(match.group("num")), CIGAROperation(match.group("code")))
            for match in self._cigar_pattern.finditer(self.cigar)
        )

    @property
    def identity(self) -> float:
        """"""
        return sum(
            num for num, code in self.iter_cigar() if code is CIGAROperation.match
        ) / min(len(self.query_sequence), len(self.target_sequence))

    def query_location(self, target_location: int) -> int:
        assert self.target_begin <= target_location < self.target_end, (
            f"The desired target location {target_location} is not within the aligned "
            f"region [{self.target_begin}; {self.target_end})."
        )
        return self.query_begin + target_location - self.target_begin
