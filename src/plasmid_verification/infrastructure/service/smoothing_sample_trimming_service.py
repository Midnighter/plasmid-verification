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


"""Provide a service that trims samples based on smoothed Phred quality."""

from typing import Tuple

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from plasmid_verification.domain.model import Sample
from plasmid_verification.domain.service import SampleTrimmingService


class SmoothingSampleTrimmingService(SampleTrimmingService):
    """Define a service that trims samples based on the smoothed Phred quality."""

    @classmethod
    def trim(
        cls,
        sample: Sample,
        *,
        prefix: str = "",
        suffix: str = "_trimmed",
        window: int = 20,
        cutoff: float = 35,
        **kwargs,
    ) -> Tuple[Sample, int, int, np.ndarray]:
        """
        Trim a sequencing sample based on the smoothed quality values and a threshold.

        Args:
            sample: A sequencing sample.
            prefix:
            suffix:
            window:
            cutoff:
            **kwargs:

        Returns:
            tuple:
                Sample: The trimmed sequencing sample.
                int: The start position of the trimmed sequence with respect to the
                    original sample.
                int: The end position of the trimmed sequence with respect to the
                    original sample.
                numpy.ndarray: The scores used by the trimming method.


        """
        index = np.arange(len(sample), dtype=int)
        smoothed = lowess(
            sample.phred_quality,
            index,
            frac=window / len(index),
            is_sorted=True,
            return_sorted=False,
        )
        smoothed[smoothed < 0.0] = 0.0
        mask = smoothed >= cutoff
        # Get min/max position that has a Phred quality higher than the threshold.
        start = index[mask][0]
        end = index[mask][-1] + 1
        return (
            Sample(
                identifier=f"{prefix}{sample.identifier}{suffix}",
                sequence=sample.sequence[start:end],
                phred_quality=sample.phred_quality[start:end],
            ),
            start,
            end,
            smoothed,
        )
