#!/usr/bin/env python3

# Copyright (C) 2024 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import argparse
import os
import sys
import requests
from tqdm import tqdm

import yaml
import pprint

def is_writeable(directory):
    try:
        with open(os.path.join(directory, 'testfile'), 'w'):
            pass
    except PermissionError:
        return False
    else:
        return True

def process_data(configuration: str, data_path: str):
    with open(configuration, 'r') as f:
        config = (yaml.safe_load(f))

        full_data_path=os.path.join(data_path, config['dataset'])

        if is_writeable(data_path):
            print("Writing data to {}".format(data_path))

            # create folder, e.g. h19
            if not os.path.exists(full_data_path):
                os.makedirs(full_data_path)

            for item in config:

                if 'url' in config[item]:
                    url = config[item]['url']
                    fname = config[item]['name']
                    type = config[item]['type']
                    with requests.get(url, stream=True) as r:
                        r.raise_for_status()
                        total_size = int(r.headers.get('content-length', 0))
                        with open(os.path.join(full_data_path, fname), 'wb') as f:
                            with tqdm(total=total_size, unit='B',
                                      unit_scale=True, desc=fname) as pbar:
                                for chunk in r.iter_content(chunk_size=8192):
                                    f.write(chunk)
                                    pbar.update(len(chunk))


def main():

    parser = argparse.ArgumentParser(
        description="Download required data from Nanopore circRNA pipeline")

    # REQUIRED ARGUMENTS
    group = parser.add_argument_group("Required options")

    group.add_argument("-c",
                       "--config",
                       dest="config",
                       help="Path to the configuration file, e.g. hg19.yml",
                       required=True
                       )

    group.add_argument("-p",
                       "--path",
                       dest="path",
                       help="Path to store the downloaded data files",
                       required=True
                       )

    args = parser.parse_args(sys.argv[1:])

    if not os.path.isfile(args.config):
        print("Configuration file not accessible.")
    else:
        process_data(configuration=args.config, data_path=args.path)

if __name__ == "__main__":
    main()
