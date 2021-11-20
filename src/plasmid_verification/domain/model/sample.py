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

from typing import Iterable, Dict, Optional

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Sample:
    def __init__(
        self,
        *,
        identifier: str,
        sequence: Seq,
        phred_quality: Iterable,
        peak_positions: Optional[Iterable] = None,
        traces: Optional[Dict[str, Iterable]] = None,
        **kwargs,
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
        if peak_positions is not None:
            self._peak_positions = np.asarray(peak_positions)
            assert len(self._sequence) == len(self._peak_positions), (
                f"The sample's ({identifier}) number of peaks must correspond to the "
                f"sequence length."
            )
        else:
            self._peak_positions = None
        if traces is not None:
            self._traces = {key: np.asarray(trace) for key, trace in traces.items()}
            assert len(self._traces) == 4, "There must be exactly one trace per base."
            for base, trace in self._traces.items():
                assert len(trace) >= len(self._peak_positions), (
                    f"Trace {base} should contain more data points than the number of "
                    f"peak positions."
                )
        else:
            self._traces = None

    @classmethod
    def from_ab1(cls, record: SeqRecord) -> Sample:
        """"""
        base_order = record.annotations["abif_raw"]["FWO_1"].decode("ascii")
        analyzed_data_channels = ("DATA9", "DATA10", "DATA11", "DATA12")
        return cls(
            identifier=record.id,
            sequence=record.seq,
            phred_quality=record.letter_annotations["phred_quality"],
            peak_positions=record.annotations["abif_raw"]["PLOC2"],
            traces={
                base: record.annotations["abif_raw"][channel]
                for base, channel in zip(base_order, analyzed_data_channels)
            },
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

    def reverse_complement(self) -> Sample:
        """Return the reverse complement of itself."""
        return Sample(
            identifier=f"complement({self.identifier})",
            sequence=self.sequence.reverse_complement(),
            phred_quality=self.phred_quality[::-1].copy(),
        )
