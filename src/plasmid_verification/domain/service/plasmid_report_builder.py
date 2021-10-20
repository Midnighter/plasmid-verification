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
)
from .sample_report_builder import SampleReportBuilder
from .sequence_alignment_service import SequenceAlignmentService


logger = logging.getLogger(__name__)


class PlasmidReportBuilder:
    def __init__(self, *, alignment_service: Type[SequenceAlignmentService], **kwargs):
        super().__init__(**kwargs)
        self.alignment_service = alignment_service
        self.report: Optional[PlasmidReport] = None

    def build(
        self,
        plasmid: Plasmid,
        samples: Iterable[Sample],
        quality_threshold: float,
        window: int,
        **kwargs,
    ) -> PlasmidReport:
        self.report = PlasmidReport(
            plasmid=plasmid,
            samples=samples,
            quality_threshold=quality_threshold,
            window=window,
        )
        if self._build_sample_reports():
            return self.report
        if self._build_summary():
            return self.report
        if self._build_effects():
            return self.report
        return self.report

    def _build_effects(self) -> bool:
        for sample in self.report.sample_reports:
            for conflict in sample.conflicts:
                conflict.evaluate_status()
                self.report.evaluate_effect(conflict)
        return False

    def _build_sample_reports(self) -> bool:
        sample_report_builder = SampleReportBuilder(
            alignment_service=self.alignment_service
        )
        for sample in self.report.samples:
            try:
                self.report.sample_reports.append(
                    sample_report_builder.build(
                        self.report.plasmid,
                        sample,
                        quality_threshold=self.report.quality_threshold,
                        window=self.report.window,
                    )
                )
            except AssertionError as error:
                logger.debug("", exc_info=error)
                msg = (
                    f"Failed to generate report for sample {sample.identifier} with "
                    f"plasmid {self.report.plasmid.identifier}."
                )
                logger.error(msg)
                self.report.errors.append(msg)
        return not self.report.sample_reports

    def _build_summary(self) -> bool:
        if len(self.report.sample_reports) < 2:
            return False
        for sample_a in self.report.sample_reports:
            if sample_a.alignment is None:
                logger.warning(
                    "Sample %s was not aligned, ignored in pairwise summary.",
                    sample_a.sample.identifier,
                )
                continue
            for sample_b in self.report.sample_reports:
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
                location = sample_b.alignment.query_location(conflict_a.begin)
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
        first_location = first_sample.alignment.query_location(first_conflict.begin)
        second_location = second_sample.alignment.query_location(second_conflict.begin)
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
        return max(first.begin, second.alignment.target_begin) <= min(
            first.end, second.alignment.target_end
        )

    @classmethod
    def _do_conflicts_overlap(cls, first: Conflict, second: Conflict) -> bool:
        return max(first.begin, second.begin) <= min(first.end, second.end)

    def _is_reliable(self, region: np.ndarray) -> bool:
        return np.nanmean(region) >= self.report.quality_threshold
