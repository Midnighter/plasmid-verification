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


"""Provide a plasmid model."""


from __future__ import annotations

from typing import List, Tuple

from .plasmid import Plasmid


class PlasmidValidationSpecification:
    @classmethod
    def is_satisfied_by(cls, plasmid: Plasmid) -> Tuple[bool, List[str], List[str]]:
        errors = []
        warnings = []
        return (
            cls._length_specification(plasmid, errors, warnings)
            and cls._feature_position_specification(plasmid, errors, warnings),
            errors,
            warnings,
        )

    @classmethod
    def _length_specification(cls, plasmid: Plasmid, errors: List[str], _) -> bool:
        if len(plasmid) == 0:
            errors.append(f"Plasmid {plasmid.identifier} has no genetic sequence.")
            return False
        return True

    @classmethod
    def _feature_position_specification(
        cls, plasmid: Plasmid, _, warnings: List[str]
    ) -> bool:
        if len(plasmid.features) == 0:
            warnings.append(f"Plasmid {plasmid.identifier} has no features.")
        for feature in plasmid.features:
            if not (
                (0 <= feature.location.start <= len(plasmid))
                and (0 <= feature.location.end <= len(plasmid))
            ):
                warnings.append(
                    f"Plasmid {plasmid.identifier} feature {feature.type} "
                    f"{feature.location.start}..{feature.location.end} is outside the "
                    f"plasmid sequence {len(plasmid)}."
                )
        return True
