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


"""Provide a sequence alignment service based on scikit-bio."""


from Bio.Seq import Seq
from skbio.alignment import StripedSmithWaterman

from plasmid_verification.domain.model import SequenceAlignment
from plasmid_verification.domain.service import SequenceAlignmentService


class ScikitBioSequenceAlignmentService(SequenceAlignmentService):
    """Define a sequence alignment service based on scikit-bio."""

    @classmethod
    def align(
        cls,
        query: Seq,
        target: Seq,
        gap_open_penalty: float = 2.0,
        gap_extension_penalty: float = 10.0,
        **kwargs,
    ) -> SequenceAlignment:
        """Return a local alignment of two given sequences."""
        align = StripedSmithWaterman(
            query_sequence=str(query),
            gap_open_penalty=int(gap_open_penalty),
            gap_extend_penalty=int(gap_extension_penalty),
        )(str(target))
        return SequenceAlignment(
            score=float(align.optimal_alignment_score),
            query_sequence=query,
            query_begin=align.query_begin,
            query_end=align.query_end,
            target_sequence=target,
            target_begin=align.target_begin,
            target_end=align.target_end_optimal,
            cigar=align.cigar,
            aligned_query_sequence=Seq(align.aligned_query_sequence),
            aligned_target_sequence=Seq(align.aligned_target_sequence),
        )
