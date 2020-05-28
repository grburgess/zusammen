from setuptools import setup
import os

import versioneer


# Create list of data files
def find_stan_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            if ".stan" in filename:
            
                paths.append(os.path.join("..", path, filename))

    return paths


stan_files = find_stan_files("zusammen/stan_models")

data_files = {"": stan_files}

setup(
    version=versioneer.get_version(),
    include_package_data=True,
    package_data=data_files,
    cmdclass=versioneer.get_cmdclass(),
)
