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


"""Provide a sequencing sample report."""


from typing import Iterable

from Bio.SeqFeature import SeqFeature

from .plasmid import Plasmid
from .sample import Sample
from .sample_report import Conflict, SampleReport


class PlasmidReport:
    def __init__(
        self,
        *,
        plasmid: Plasmid,
        samples: Iterable[Sample] = (),
        quality_threshold: float,
        sample_reports: Iterable[SampleReport] = (),
        errors: Iterable[str] = (),
        warnings: Iterable[str] = (),
        **kwargs
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self.plasmid = plasmid
        self.samples = list(samples)
        self.quality_threshold = quality_threshold
        self.sample_reports = list(sample_reports)
        self.errors = list(errors)
        self.warnings = list(warnings)

    def evaluate_effect(self, conflict: Conflict) -> None:
        for feature in self.plasmid.features:
            if feature.type != "CDS":
                continue
            if not self.do_feature_and_conflict_overlap(feature, conflict):
                continue
            # TODO: Migrate code from sanger-sequencing.

    @classmethod
    def do_feature_and_conflict_overlap(
        cls, feature: SeqFeature, conflict: Conflict
    ) -> bool:
        return max(feature.location.start, conflict.plasmid_begin) <= min(
            feature.location.plasmid_end, conflict.plasmid_end
        )
