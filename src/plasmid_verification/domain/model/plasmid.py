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

import logging
from typing import Dict, Iterable, List, Type

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .sample import Sample
from .sequence_alignment import CIGAROperation
from ..service import SampleTrimmingService, SequenceAlignmentService


logger = logging.getLogger(__name__)


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

    def verify(
        self,
        sample: Sample,
        trimming_service: Type[SampleTrimmingService],
        alignment_service: Type[SequenceAlignmentService],
        **kwargs,
    ) -> None:
        """"""
        trimmed = trimming_service.trim(sample, **kwargs)
        alignment = alignment_service.align(
            query=trimmed.sequence, target=self.sequence, **kwargs
        )
        # reverse = alignment_service.align(
        #     query=trimmed.sequence.reverse_complement(),
        #     target=self.sequence,
        #     **kwargs,
        # )
        if False:
            # if reverse.identity > alignment.identity:
            alignment = reverse
            strand = -1
            note = (
                f"Alignment of reverse complement of trimmed sample "
                f"{sample.identifier} to plasmid {self.identifier}."
            )
        else:
            strand = 1
            note = (
                f"Alignment of trimmed sample {sample.identifier} to plasmid "
                f"{self.identifier}."
            )
        logger.info(
            "Add sample from %d to %d.", alignment.target_begin, alignment.target_end
        )
        self.add_feature(
            alignment.target_begin,
            alignment.target_end,
            type="misc_feature",
            strand=strand,
            qualifiers={
                "label": f"Verified by {sample.identifier}",
                "note": note,
            },
        )
        index = alignment.target_begin
        for num, code in alignment.iter_cigar():
            if code is CIGAROperation.match:
                index += num
                continue
            logger.info("Add conflict from %d to %d.", index, index + num)
            self.add_feature(
                index,
                index + num,
                type="conflict",
                strand=strand,
                qualifiers={
                    "label": str(code.value),
                },
            )
            index += num
