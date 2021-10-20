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


""""""


import logging
from typing import Iterable, Type

from plasmid_verification.domain.model import (
    PlasmidReport,
    PlasmidRepository,
    SampleRepository,
)
from plasmid_verification.domain.service import (
    PlasmidReportBuilder,
    SequenceAlignmentService,
)


logger = logging.getLogger(__name__)


class PlasmidReportService:
    def __init__(
        self,
        *,
        plasmid_repository: PlasmidRepository,
        sample_repository: SampleRepository,
        alignment_service: Type[SequenceAlignmentService],
        **kwargs,
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self.plasmid_repository = plasmid_repository
        self.sample_repository = sample_repository
        self._builder = PlasmidReportBuilder(alignment_service=alignment_service)

    def report(
        self,
        plasmid_id: str,
        sample_ids: Iterable[str],
        quality_threshold: float,
        window: int,
    ) -> PlasmidReport:
        plasmid = self.plasmid_repository.find(plasmid_id)
        assert plasmid is not None, f"Plasmid {plasmid_id} not found in repository."
        samples = []
        errors = []
        for sample_id in sample_ids:
            sample = self.sample_repository.find(sample_id)
            if sample is None:
                errors.append(f"Sample {sample_id} not found in repository.")
            else:
                samples.append(sample)
        result = self._builder.build(
            plasmid=plasmid,
            samples=samples,
            quality_threshold=quality_threshold,
            window=window,
        )
        result.errors.extend(errors)
        return result
