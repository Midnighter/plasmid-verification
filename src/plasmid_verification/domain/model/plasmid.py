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


"""Provide a plasmid model."""


from __future__ import annotations

from typing import Dict, Iterable, List

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


class Plasmid:
    def __init__(
        self,
        *,
        identifier: str,
        sequence: Seq,
        features: Iterable[SeqFeature],
        annotation: dict,
        **kwargs,
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self._identifier = identifier
        self._sequence = sequence
        self._features = sorted(
            features, key=lambda f: (f.location.start, f.location.end)
        )
        self._annotation = annotation

    @classmethod
    def from_genbank(cls, record: SeqRecord) -> Plasmid:
        """"""
        return cls(
            identifier=record.id,
            sequence=record.seq,
            features=record.features,
            annotation=record.annotations,
        )

    def __len__(self) -> int:
        """"""
        return len(self.sequence)

    @property
    def identifier(self) -> str:
        """Return the plasmid's identifier."""
        return self._identifier

    @property
    def sequence(self) -> Seq:
        """Return the DNA sequence of the plasmid."""
        return self._sequence

    @property
    def features(self) -> List[SeqFeature]:
        """Return the features annotated on this plasmid."""
        return self._features.copy()

    @property
    def annotation(self) -> dict:
        """Return the plasmid annotation."""
        return self._annotation.copy()

    def add_feature(
        self,
        start: int,
        end: int,
        qualifiers: Dict[str, str],
        type: str,
        strand: int = 1,
    ) -> None:
        self._features.append(
            SeqFeature(
                location=FeatureLocation(start=start, end=end, strand=strand),
                type=type,
                qualifiers=qualifiers,
            )
        )
        self._features.sort(key=lambda f: (f.location.start, f.location.end))

    def to_seqrecord(self) -> SeqRecord:
        return SeqRecord(
            id=self.identifier,
            seq=self.sequence,
            features=self.features,
            annotations=self.annotation,
        )
