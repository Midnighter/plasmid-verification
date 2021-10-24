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


"""Provide an abstract interface for a sample trimming service."""


from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np

from ..model import Sample


class SampleTrimmingService(ABC):
    """Define the abstract interface for a sample trimming service."""

    @classmethod
    @abstractmethod
    def trim(
        cls,
        sample: Sample,
        *,
        cutoff: float,
        prefix: str = "",
        suffix: str = "_trimmed",
        **kwargs
    ) -> Tuple[Sample, int, int, np.ndarray]:
        """
        Trim a sequencing sample.

        Args:
            sample: A sequencing sample.
            cutoff:
            prefix:
            suffix:
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
