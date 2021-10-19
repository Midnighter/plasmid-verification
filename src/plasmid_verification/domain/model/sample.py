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


"""Provide a sequencing sample model."""


from __future__ import annotations

from typing import Iterable

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Sample:
    def __init__(
        self, *, identifier: str, sequence: Seq, phred_quality: Iterable, **kwargs
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self._identifier = identifier
        self._sequence = sequence
        self._phred_quality = np.asarray(phred_quality)
        assert len(self._sequence) == len(self._phred_quality), (
            f"The sample's ({identifier}) sequence must have a Phred quality score for "
            f"each position."
        )

    @classmethod
    def from_ab1(cls, record: SeqRecord) -> Sample:
        """"""
        return cls(
            identifier=record.id,
            sequence=record.seq,
            phred_quality=record.letter_annotations["phred_quality"],
        )

    def __len__(self) -> int:
        """"""
        return len(self.sequence)

    @property
    def identifier(self) -> str:
        """Return the sample's identifier."""
        return self._identifier

    @property
    def sequence(self) -> Seq:
        """Return the DNA sequence of the sample."""
        return self._sequence

    @property
    def phred_quality(self) -> np.ndarray:
        """Return the Phred quality scores for all sequence positions."""
        return self._phred_quality.copy()
