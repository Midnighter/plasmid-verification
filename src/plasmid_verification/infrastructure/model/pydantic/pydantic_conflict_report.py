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

from pydantic import Field, PositiveInt, conint

from plasmid_verification.domain.model import Conflict
from .pydantic_custom_base import PydanticCustomBase
from plasmid_verification.domain.model import (
    ConflictType,
    ConflictReliability,
    ConflictStatus,
)


class PydanticConflictReport(PydanticCustomBase):

    plasmid_letters: str = Field(..., alias="plasmidLetters", title="Plasmid Letters")
    sample_letters: str = Field(..., alias="sampleLetters", title="Sample Letters")
    plasmid_position: PositiveInt = Field(
        ..., alias="plasmidPosition", title="Plasmid Position"
    )
    span: PositiveInt
    # sample_position: PositiveInt = Field(
    #     ..., alias="samplePosition", title="Sample Position"
    # )
    type: ConflictType
    reliability: ConflictReliability
    num_confirmed: conint(ge=0) = Field(
        ..., alias="numConfirmed", title="Number Confirmed"
    )
    num_invalidated: conint(ge=0) = Field(
        ..., alias="numInvalidated", title="Number Invalidated"
    )
    status: ConflictStatus
    # features_hit: List[SequenceFeature] = Field(
    #     (), alias="featuresHit", title="Hit Features"
    # )
    # effects: List[Effect] = ()

    class Config(PydanticCustomBase.Config):

        description = "Summarize a detected alignment conflict."

    @classmethod
    def from_domain(cls, conflict: Conflict) -> PydanticConflictReport:
        result = cls(
            plasmid_letters="",
            sample_letters="",
            plasmid_position=conflict.begin,
            span=conflict.span,
            type=conflict.type,
            reliability=conflict.reliability,
            num_confirmed=conflict.num_confirmed,
            num_invalidated=conflict.num_invalidated,
            status=conflict.status,
        )
        return result
