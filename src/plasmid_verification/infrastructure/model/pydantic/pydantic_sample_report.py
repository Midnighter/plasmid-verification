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


from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional

from pydantic import Field, PositiveInt, confloat, conint

from plasmid_verification.domain.model import ConflictStatus, SampleReport

from .pydantic_conflict_report import PydanticConflictReport
from .pydantic_custom_base import PydanticCustomBase


class PydanticSampleReport(PydanticCustomBase):

    id: str = Field(..., description="The given identifier for the sample.")
    read_length: PositiveInt = Field(
        ...,
        alias="readLength",
        title="Read Length",
        description="The number of nucleotides of this sample read.",
    )
    median_quality: confloat(ge=0.0, le=65.0) = Field(
        ...,
        alias="medianQuality",
        title="Median Quality",
        description="The median Phred quality of the sample read.",
    )
    trim_start: Optional[conint(ge=0)] = Field(
        None,
        alias="trimStart",
        title="Trim Start",
        description="The number of nucleotides trimmed at the beginning of the "
        "sequence due to low Phred quality and before alignment.",
    )
    trim_end: Optional[conint(ge=0)] = Field(
        None,
        alias="trimEnd",
        title="Trim End",
        description="The number of nucleotides trimmed at the end of the sequence due "
        "to low Phred quality and before alignment.",
    )
    phred_quality: List[confloat(ge=0.0, le=65.0)] = Field(
        ...,
        alias="phredQuality",
        title="Phred Quality",
        description="The Phred quality score at each position of the sequence.",
    )
    smoothed_scores: Optional[List[confloat(ge=0.0, le=65.0)]] = Field(
        None,
        alias="smoothedScores",
        title="Smoothed Scores",
        description="The smoothed Phred quality scores.",
    )
    sequence: str = Field(..., description="The full nucleotide sequence.")
    strand: Optional[int] = None
    conflicts: List[PydanticConflictReport] = Field(
        (), description="A collection of individual alignment conflict reports."
    )
    conflict_counts: Dict[ConflictStatus, int] = Field(
        {}, alias="conflictCounts", title="Conflict Counts"
    )
    errors: List[str] = Field(
        (), description="Any errors that occurred during sample analysis."
    )
    warnings: List[str] = Field(
        (), description="Any warnings that occurred during sample analysis."
    )

    class Config(PydanticCustomBase.Config):

        description = "Summarize the results for a sample sequencing read."

    @classmethod
    def from_domain(cls, report: SampleReport) -> PydanticSampleReport:
        result = cls(
            id=report.sample.identifier,
            read_length=len(report.sample),
            median_quality=report.median_quality,
            trim_start=report.trim_begin,
            trim_end=report.trim_end,
            phred_quality=list(report.sample.phred_quality),
            sequence=str(report.sample.sequence),
            strand=report.strand,
            errors=report.errors,
            warnings=report.warnings,
        )
        if report.alignment is not None:
            result.smoothed_scores = list(report.smoothed)
        counter = Counter()
        for conflict in report.conflicts:
            result.conflicts.append(PydanticConflictReport.from_domain(conflict))
            counter[conflict.status] += 1
        result.conflict_counts = {status: counter[status] for status in ConflictStatus}
        return result
