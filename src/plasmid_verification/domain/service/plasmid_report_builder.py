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


import logging
from typing import Iterable, Optional, Type

import numpy as np

from ..model import (
    Conflict,
    ConflictReliability,
    Plasmid,
    PlasmidReport,
    Sample,
    SampleReport,
    ConflictStatus,
    StrandDirection,
    ConflictType,
)
from .sample_report_builder import SampleReportBuilder
from .sample_trimming_service import SampleTrimmingService
from .sequence_alignment_service import SequenceAlignmentService


logger = logging.getLogger(__name__)


class PlasmidReportBuilder:
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
        self._report: Optional[PlasmidReport] = None

    def build(
        self,
        plasmid: Plasmid,
        samples: Iterable[Sample],
        quality_threshold: float,
        trim_kwargs: dict,
        alignment_kwargs: dict,
        **kwargs,
    ) -> PlasmidReport:
        self._report = PlasmidReport(
            plasmid=plasmid,
            samples=samples,
            quality_threshold=quality_threshold,
            **kwargs,
        )
        if self._build_sample_reports(trim_kwargs, alignment_kwargs):
            return self._report
        if self._build_summary():
            return self._report
        if self._build_effects():
            return self._report
        return self._report

    def _build_sample_reports(
        self,
        trim_kwargs: dict,
        alignment_kwargs: dict,
    ) -> bool:
        sample_report_builder = SampleReportBuilder(
            trimming_service=self._trimming_service,
            alignment_service=self._alignment_service,
        )
        for sample in self._report.samples:
            try:
                self._report.sample_reports.append(
                    sample_report_builder.build(
                        self._report.plasmid,
                        sample,
                        quality_threshold=self._report.quality_threshold,
                        trim_kwargs=trim_kwargs,
                        alignment_kwargs=alignment_kwargs,
                    )
                )
            except AssertionError as error:
                logger.debug("", exc_info=error)
                msg = (
                    f"Failed to generate report for sample {sample.identifier} with "
                    f"plasmid {self._report.plasmid.identifier}."
                )
                logger.error(msg)
                self._report.errors.append(msg)
        return not self._report.sample_reports

    def _build_summary(self) -> bool:
        if len(self._report.sample_reports) < 2:
            return False
        for sample_a in self._report.sample_reports:
            if sample_a.alignment is None:
                logger.warning(
                    "Sample %s was not aligned, ignored in pairwise summary.",
                    sample_a.sample.identifier,
                )
                continue
            for sample_b in self._report.sample_reports:
                if sample_a is sample_b:
                    continue
                if sample_b.alignment is None:
                    logger.warning(
                        "Sample %s was not aligned, ignored in pairwise summary.",
                        sample_a.sample.identifier,
                    )
                    continue
                if not self._do_samples_overlap(sample_a, sample_b):
                    logger.debug(
                        "Sample %s and %s do not overlap.",
                        sample_a.sample.identifier,
                        sample_b.sample.identifier,
                    )
                    continue
                self._compare_conflicts(sample_a, sample_b)
        return False

    def _compare_conflicts(self, sample_a: SampleReport, sample_b: SampleReport):
        for conflict_a in sample_a.conflicts:
            if not self._do_sample_and_conflict_overlap(conflict_a, sample_b):
                continue
            for conflict_b in sample_b.conflicts:
                if self._do_conflicts_overlap(conflict_a, conflict_b):
                    if conflict_b.reliability == ConflictReliability.LOW:
                        break
                    if conflict_a.type != conflict_b.type:
                        conflict_a.invalidate()
                    elif self._equal_sequence(
                        conflict_a, sample_a, conflict_b, sample_b
                    ):
                        conflict_a.confirm()
                    else:
                        conflict_a.invalidate()
                    break
            else:
                # Can sample B invalidate conflict A?
                location = sample_b.alignment.query_location(conflict_a.plasmid_begin)
                if self._is_reliable(
                    sample_b.trimmed.phred_quality[
                        location : location + conflict_a.span
                    ]
                ):
                    conflict_a.invalidate()

    @classmethod
    def _equal_sequence(
        cls,
        first_conflict: Conflict,
        first_sample: SampleReport,
        second_conflict: Conflict,
        second_sample: SampleReport,
    ) -> bool:
        first_location = first_sample.alignment.query_location(
            first_conflict.plasmid_begin
        )
        second_location = second_sample.alignment.query_location(
            second_conflict.plasmid_begin
        )
        return (
            first_sample.trimmed.sequence[
                first_location : first_location + first_conflict.span
            ]
            == second_sample.trimmed.sequence[
                second_location : second_location + second_conflict.span
            ]
        )

    @classmethod
    def _do_samples_overlap(cls, first: SampleReport, second: SampleReport) -> bool:
        return max(first.alignment.target_begin, second.alignment.target_begin) <= min(
            first.alignment.target_end, second.alignment.target_end
        )

    @classmethod
    def _do_sample_and_conflict_overlap(
        cls, first: Conflict, second: SampleReport
    ) -> bool:
        return max(first.plasmid_begin, second.alignment.target_begin) <= min(
            first.plasmid_end, second.alignment.target_end
        )

    @classmethod
    def _do_conflicts_overlap(cls, first: Conflict, second: Conflict) -> bool:
        return max(first.plasmid_begin, second.plasmid_begin) <= min(
            first.plasmid_end, second.plasmid_end
        )

    def _is_reliable(self, region: np.ndarray) -> bool:
        return np.nanmean(region) >= self._report.quality_threshold

    def _build_effects(self) -> bool:
        for sample in self._report.sample_reports:
            for conflict in sample.conflicts:
                conflict.evaluate_status()
                self._report.evaluate_effect(conflict)
                if (
                    conflict.status is not ConflictStatus.RESOLVED
                    and conflict.reliability is ConflictReliability.HIGH
                ):
                    self._report.plasmid.add_feature(
                        conflict.plasmid_begin,
                        conflict.plasmid_end,
                        type="conflict",
                        strand=StrandDirection.FORWARD,
                        qualifiers={
                            "label": self._conflict_label(conflict, sample),
                            "note": "Effect: #TODO",
                        },
                    )
        return False

    def _conflict_label(self, conflict: Conflict, sample: SampleReport) -> str:
        end = (
            conflict.plasmid_end + 1
            if conflict.plasmid_begin == conflict.plasmid_end
            else conflict.plasmid_end
        )
        plasmid_seq = self._report.plasmid.sequence[conflict.plasmid_begin : end]
        end = (
            conflict.sample_end + 1
            if conflict.sample_begin == conflict.sample_end
            else conflict.sample_end
        )
        sample_seq = sample.sample.sequence[conflict.sample_begin : end]
        if conflict.type is ConflictType.VARIANT:
            label = f"SNV: {plasmid_seq} -> {sample_seq}"
        elif conflict.type is ConflictType.DELETION:
            label = f"Deletion: {plasmid_seq} -> {'-' * conflict.span}"
        elif conflict.type is ConflictType.INSERTION:
            label = f"Deletion: {'-' * conflict.span} -> {sample_seq}"
        else:
            raise RuntimeError(f"Unknown conflict type {conflict.type}.")
        return label
