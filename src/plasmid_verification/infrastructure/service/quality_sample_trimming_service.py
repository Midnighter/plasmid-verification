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


"""Provide a service that trims samples based on Phred quality."""


import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from plasmid_verification.domain.model import Sample
from plasmid_verification.domain.service import SampleTrimmingService


class QualitySampleTrimmingService(SampleTrimmingService):
    """Define a service that trims samples based on Phred quality."""

    @classmethod
    def trim(
        cls,
        sample: Sample,
        quality_threshold: float = 35,
        window: int = 20,
        prefix: str = "",
        suffix: str = "_trimmed",
        **kwargs,
    ) -> Sample:
        """Trim a sequence sample based on quality."""
        index = np.arange(len(sample), dtype=int)
        smoothed = lowess(
            sample.phred_quality, index, is_sorted=True, return_sorted=False
        )
        mask = smoothed >= quality_threshold
        min_i = index[mask][0]
        max_i = index[mask][-1] + 1
        return Sample(
            identifier=f"{prefix}{sample.identifier}{suffix}",
            sequence=sample.sequence[min_i:max_i],
            phred_quality=sample.phred_quality[min_i:max_i],
        )
