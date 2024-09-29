from setuptools import setup, find_packages
import re

# Read the contents of README.md
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read the contents of fabbit/main.py to extract the version
with open("fabbit/main.py", "r", encoding="utf-8") as f:
    content = f.read()
    version_match = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", content, re.M)
    if version_match:
        version = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string in fabbit/main.py")

setup(
    name="fabbit",
    version=version,  # Use the version extracted from fabbit/main.py
    author="Your Actual Name",
    author_email="your.actual@email.com",
    description="Fabbit: FAst coregenome alignment Based on Bidirectional best hIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EnzoAndree/FABBIT",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        # Add your dependencies here, for example:
        'biopython',
        'numpy',
        'pandas',
        'matplotlib',
        'tqdm',
        'pyrodigal',
        'pyfastx',
    ],
    entry_points={
        "console_scripts": [
            "fabbit=fabbit.main:main",
        ],
    },
)