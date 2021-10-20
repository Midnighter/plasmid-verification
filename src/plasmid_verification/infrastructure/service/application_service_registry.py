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

from plasmid_verification.application.service import SampleReportService
from .domain_registry import DomainRegistry


class ApplicationServiceRegistry:
    def __init__(self, domain_registry: Type[DomainRegistry], **kwargs) -> None:
        super().__init__(**kwargs)
        self.domain_registry = domain_registry

    def sample_report_service(
        self, archive: Optional[zipfile.ZipFile] = None
    ) -> SampleReportService:
        return SampleReportService(
            plasmid_repository=self.domain_registry.plasmid_repository(archive=archive),
            sample_repository=self.domain_registry.sample_repository(archive=archive),
            trimming_service=self.domain_registry.sample_trimming_service(),
            alignment_service=self.domain_registry.sequence_alignment_service(),
        )
