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

from io import StringIO
from typing import List

from Bio import SeqIO
from pydantic import Field

from plasmid_verification.domain.model import PlasmidReport
from .pydantic_custom_base import PydanticCustomBase
from .pydantic_sample_report import PydanticSampleReport


class PydanticPlasmidReport(PydanticCustomBase):

    id: str = Field(..., description="The given identifier for the plasmid.")
    genbank: str = Field(
        ..., description="The plasmid and all its feature in Genbank format."
    )
    samples: List[PydanticSampleReport] = Field(
        (), description="A collection of individual sample (sequencing read) reports."
    )
    errors: List[str] = Field(
        (), description="Any errors that occurred during plasmid analysis."
    )
    warnings: List[str] = Field(
        (), description="Any warnings that occurred during plasmid analysis."
    )

    class Config(PydanticCustomBase.Config):

        description = "Summarize the results for an entire plasmid."

    @classmethod
    def from_domain(cls, report: PlasmidReport) -> PydanticPlasmidReport:
        genbank = StringIO()
        SeqIO.write(report.plasmid.to_seqrecord(), genbank, "genbank")
        result = cls(
            id=report.plasmid.identifier,
            genbank=genbank.getvalue(),
            errors=report.errors,
            warnings=report.warnings,
        )
        for sample in report.sample_reports:
            result.samples.append(PydanticSampleReport.from_domain(sample))
        return result
