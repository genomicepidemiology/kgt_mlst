import yaml
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '')] + sys.path

import kgt_mlst.version as version

data = {
    "package": {
        "name": "kgt_mlst",
        "version": version.__version__
    },
    "source": {
        "url": "https://github.com/genomicepidemiology/kgt_mlst/archive/refs/tags/{}.tar.gz".format(version.__version__),
    },
    "build": {
        "number": 0,
        "noarch": "python",
        "script": "{{ PYTHON }} -m pip install . --no-deps --ignore-installed -vvv"
    },
    "requirements": {
        "host": [
            "python >=3.6",
            "pip"
        ],
        "run": [
            "python >=3.6",
            "kma >=1.4.9"
        ]
    },
    "about": {
        "home": "https://github.com/genomicepidemiology/kgt_mlst",
        "summary": "kgt_mlst test.",
        "license": "Apache-2.0"
    }
}

# Convert the data to YAML and print it
os.system('mkdir conda')
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
