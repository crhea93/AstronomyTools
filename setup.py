import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AstronomyTools",
    version="0.0.1",
    author="Carter Rhea",
    author_email="carterrhea93@gmail.com",
    description="Collection of Useful Astronomy Tools that I wrote",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crhea93/AstronomyTools",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: General Use License",
        "Operating System :: Linux OS",
    ),
)
