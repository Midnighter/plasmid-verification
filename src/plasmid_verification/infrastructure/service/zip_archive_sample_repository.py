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


"""Provide a sample repository from a zip archive."""


import zipfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from plasmid_verification.domain.model import Sample, SampleRepository


class ZipArchiveSampleRepository(SampleRepository):
    """Define a sample repository from a zip archive."""

    _VALID_EXTENSIONS = frozenset((".ab1", ".abi"))

    def __init__(self, *, archive: zipfile.ZipFile, **kwargs):
        """"""
        super().__init__(**kwargs)
        self._samples = {
            path.stem: Sample.from_ab1(
                SeqIO.read(archive.open(str(path)), format="abi")
            )
            for path in (Path(filename) for filename in archive.namelist())
            if path.suffix in self._VALID_EXTENSIONS
        }

    def add(self, sample: Sample) -> None:
        """"""
        self._samples[sample.identifier] = sample

    def find(self, identifier: str) -> Optional[Sample]:
        """"""
        return self._samples.get(identifier)
