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
from typing import Type

import numpy as np

from plasmid_verification.domain.model import (
    Plasmid,
    Sample,
    PlasmidRepository,
    SampleRepository,
    CIGAROperation,
)
from plasmid_verification.domain.service import (
    SampleTrimmingService,
    SequenceAlignmentService,
)


logger = logging.getLogger(__name__)


class SampleReportService:
    def __init__(
        self,
        *,
        plasmid_repository: PlasmidRepository,
        sample_repository: SampleRepository,
        trimming_service: Type[SampleTrimmingService],
        alignment_service: Type[SequenceAlignmentService],
        **kwargs,
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self.plasmid_repository = plasmid_repository
        self.sample_repository = sample_repository
        self.trimming_service = trimming_service
        self.alignment_service = alignment_service

    def report(
        self,
        plasmid_id: str,
        sample_id: str,
        quality_threshold: float = 35,
        window: int = 20,
    ):
        plasmid = self.plasmid_repository.find(plasmid_id)
        sample = self.sample_repository.find(sample_id)
        self.verify(
            plasmid,
            sample,
            quality_threshold=quality_threshold,
            window=window,
        )

    def verify(
        self,
        plasmid: Plasmid,
        sample: Sample,
        quality_threshold: float = 35,
        window: int = 20,
        **kwargs,
    ) -> None:
        """"""
        if (median := np.nanmedian(sample.phred_quality)) < quality_threshold:
            raise AssertionError(
                f"The sample's ({sample.identifier}) median Phred quality {median} is "
                f"below the threshold {quality_threshold}."
            )
        trimmed = self.trimming_service.trim(
            sample, quality_threshold=quality_threshold, window=window, **kwargs
        )
        alignment = self.alignment_service.align(
            query=trimmed.sequence, target=plasmid.sequence, **kwargs
        )
        # reverse = alignment_service.align(
        #     query=trimmed.sequence.reverse_complement(),
        #     target=plasmid.sequence,
        #     **kwargs,
        # )
        if False:
            # if reverse.identity > alignment.identity:
            alignment = reverse
            strand = -1
            note = (
                f"Alignment of reverse complement of trimmed sample "
                f"{sample.identifier} to plasmid {plasmid.identifier}."
            )
        else:
            strand = 1
            note = (
                f"Alignment of trimmed sample {sample.identifier} to plasmid "
                f"{plasmid.identifier}."
            )
        logger.info(
            "Add sample from %d to %d.", alignment.target_begin, alignment.target_end
        )
        plasmid.add_feature(
            alignment.target_begin,
            alignment.target_end,
            type="misc_feature",
            strand=strand,
            qualifiers={
                "label": f"Verified by {sample.identifier}",
                "note": note,
            },
        )
        index = alignment.target_begin
        for num, code in alignment.iter_cigar():
            if code is CIGAROperation.match:
                index += num
                continue
            logger.info("Add conflict from %d to %d.", index, index + num)
            plasmid.add_feature(
                index,
                index + num,
                type="conflict",
                strand=strand,
                qualifiers={
                    "label": str(code.value),
                },
            )
            index += num
