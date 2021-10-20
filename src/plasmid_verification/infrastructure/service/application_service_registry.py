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
from typing import Optional

from plasmid_verification.application.service import PlasmidReportService

from .domain_registry import DomainRegistry


class ApplicationServiceRegistry:
    @classmethod
    def plasmid_report_service(
        cls, archive: Optional[zipfile.ZipFile] = None
    ) -> PlasmidReportService:
        return PlasmidReportService(
            plasmid_repository=DomainRegistry.plasmid_repository(archive=archive),
            sample_repository=DomainRegistry.sample_repository(archive=archive),
            alignment_service=DomainRegistry.sequence_alignment_service(),
        )
