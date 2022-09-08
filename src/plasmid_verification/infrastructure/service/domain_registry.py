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


import zipfile
from typing import Optional, Type

from plasmid_verification.domain.model import PlasmidRepository, SampleRepository
from plasmid_verification.domain.service import (
    SampleTrimmingService,
    SequenceAlignmentService,
)

from .zip_archive_plasmid_repository import ZipArchivePlasmidRepository
from .zip_archive_sample_repository import ZipArchiveSampleRepository


class DomainRegistry:
    @classmethod
    def plasmid_repository(
        cls, archive: Optional[zipfile.ZipFile] = None
    ) -> PlasmidRepository:
        if archive is None:
            raise RuntimeError(
                "Currently, no other plasmid repository format is supported."
            )
        return ZipArchivePlasmidRepository(archive=archive)

    @classmethod
    def sample_repository(
        cls, archive: Optional[zipfile.ZipFile] = None
    ) -> SampleRepository:
        if archive is None:
            raise RuntimeError(
                "Currently, no other sample repository format is supported."
            )
        return ZipArchiveSampleRepository(archive=archive)

    @classmethod
    def sample_trimming_service(
        cls, use_smooth: bool = False
    ) -> Type[SampleTrimmingService]:
        if use_smooth:
            from .smoothing_sample_trimming_service import (
                SmoothingSampleTrimmingService,
            )

            return SmoothingSampleTrimmingService
        else:
            from .error_probability_sample_trimming_service import (
                ErrorProbabilitySampleTrimmingService,
            )

            return ErrorProbabilitySampleTrimmingService

    @classmethod
    def sequence_alignment_service(cls) -> Type[SequenceAlignmentService]:
        from plasmid_verification.infrastructure.service.skbio import (
            ScikitBioSequenceAlignmentService,
        )

        return ScikitBioSequenceAlignmentService
