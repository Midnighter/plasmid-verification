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

import logging
from typing import Optional, Type

import numpy as np

from ..model import (
    CIGAROperation,
    Conflict,
    ConflictType,
    Plasmid,
    Sample,
    SampleReport,
    StrandDirection,
)
from .sample_trimming_service import SampleTrimmingService
from .sequence_alignment_service import SequenceAlignmentService

logger = logging.getLogger(__name__)


class SampleReportBuilder:
    def __init__(
        self,
        *,
        trimming_service: Type[SampleTrimmingService],
        alignment_service: Type[SequenceAlignmentService],
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._trimming_service = trimming_service
        self._alignment_service = alignment_service
        self._report: Optional[SampleReport] = None

    def build(
        self,
        plasmid: Plasmid,
        sample: Sample,
        quality_threshold: float,
        trim_kwargs: dict,
        alignment_kwargs: dict,
        **kwargs,
    ) -> SampleReport:
        self._report = SampleReport(
            sample=sample,
            quality_threshold=quality_threshold,
            trim_kwargs=trim_kwargs,
            alignment_kwargs=alignment_kwargs,
            **kwargs,
        )
        if self._build_median():
            return self._report
        if self._build_trimmed(**trim_kwargs):
            return self._report
        if self._build_alignment(plasmid, **alignment_kwargs):
            return self._report
        if self._build_conflicts():
            return self._report
        return self._report

    def _build_median(self) -> bool:
        self._report.median_quality = np.nanmedian(self._report.sample.phred_quality)
        if self._report.median_quality < self._report.quality_threshold:
            self._report.errors.append(
                f"The sample's ({self._report.sample.identifier}) median Phred quality "
                f"{self._report.median_quality} is below the threshold "
                f"{self._report.quality_threshold}."
            )
            return True
        return False

    def _build_trimmed(
        self, prefix: str = "", suffix: str = "_trimmed", **kwargs
    ) -> bool:
        """Trim a sequence sample based on Phred quality."""
        (
            self._report.trimmed,
            self._report.trim_begin,
            self._report.trim_end,
            self._report.smoothed,
        ) = self._trimming_service.trim(
            self._report.sample,
            prefix=prefix,
            suffix=suffix,
            **kwargs,
        )
        return False

    def _build_alignment(self, plasmid: Plasmid, **kwargs) -> bool:
        self._report.alignment = self._alignment_service.align(
            query=self._report.trimmed.sequence, target=plasmid.sequence, **kwargs
        )
        self._report.strand = StrandDirection.FORWARD
        direction = "forward"
        note = (
            f"Alignment of trimmed sample {self._report.sample.identifier} to "
            f"plasmid {plasmid.identifier}."
        )
        if self._report.alignment.identity < 0.95:
            reverse_alignment = self._alignment_service.align(
                query=self._report.trimmed.sequence.reverse_complement(),
                target=plasmid.sequence,
                **kwargs,
            )
            if reverse_alignment.score > self._report.alignment.score:
                self._report.alignment = reverse_alignment
                self._report.trimmed = self._report.trimmed.reverse_complement()
                self._report.sample = self._report.sample.reverse_complement()
                self._report.strand = StrandDirection.REVERSE
                direction = "reverse"
                note = (
                    f"Alignment of reverse complement of trimmed sample "
                    f"{self._report.sample.identifier} to plasmid "
                    f"{plasmid.identifier}."
                )
        logger.info(
            "Add sample from %d to %d in %s direction.",
            self._report.alignment.target_begin,
            self._report.alignment.target_end,
            direction,
        )
        plasmid.add_feature(
            self._report.alignment.target_begin,
            self._report.alignment.target_end,
            type="misc_feature",
            strand=self._report.strand,
            qualifiers={
                "label": f"Alignment of {self._report.sample.identifier}",
                "note": note,
            },
        )
        return False

    def _build_conflicts(self) -> bool:
        target_index = self._report.alignment.target_begin
        query_index = self._report.alignment.query_begin
        events = list(self._report.alignment.iter_cigar())
        for idx, (span, code) in enumerate(events):
            if code is CIGAROperation.match:
                for q, t in zip(
                    range(query_index, query_index + span),
                    range(target_index, target_index + span),
                ):
                    if (
                        self._report.alignment.query_sequence[q]
                        != self._report.alignment.target_sequence[t]
                    ):
                        self._report.conflicts.append(
                            Conflict(
                                type=ConflictType.VARIANT,
                                plasmid_begin=t,
                                plasmid_end=t,
                                sample_begin=q,
                                sample_end=q,
                                reliability=self._report.judge_reliability(
                                    self._report.trimmed.phred_quality[q]
                                ),
                                span=1,
                            )
                        )
                target_index += span
                query_index += span
            elif code is CIGAROperation.deletion:
                # Our logic requires that the preceding and following event are matches.
                assert events[idx - 1][1] is CIGAROperation.match
                assert events[idx + 1][1] is CIGAROperation.match
                quality = (
                    self._report.trimmed.phred_quality[query_index - 1]
                    + self._report.trimmed.phred_quality[query_index]
                ) / 2
                logger.info(
                    "Add deletion from %d to %d.", target_index, target_index + span
                )
                self._report.conflicts.append(
                    Conflict(
                        type=ConflictType.DELETION,
                        plasmid_begin=target_index,
                        plasmid_end=target_index + span,
                        sample_begin=query_index,
                        sample_end=query_index,
                        reliability=self._report.judge_reliability(quality),
                        span=span,
                    )
                )
                target_index += span
            elif code is CIGAROperation.insertion:
                # Our logic requires that the preceding and following event are matches.
                assert events[idx - 1][1] is CIGAROperation.match
                assert events[idx + 1][1] is CIGAROperation.match
                quality = np.nanmean(
                    self._report.trimmed.phred_quality[query_index : query_index + span]
                )
                logger.info("Add insertion at %d.", target_index)
                self._report.conflicts.append(
                    Conflict(
                        type=ConflictType.INSERTION,
                        plasmid_begin=target_index,
                        plasmid_end=target_index,
                        sample_begin=query_index,
                        sample_end=query_index + span,
                        reliability=self._report.judge_reliability(quality),
                        span=span,
                    )
                )
                query_index += span
            else:
                raise ValueError("Unhandled CIGAR event.")
        return False
