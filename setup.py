import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='ebspy',
    version="0.0.1",
    author="Tom McAuliffe",
    author_email="tmcauliffe@icloud.com",
    description="A package for handling EBSPs",
    long_description=long_description,
    long_description_content_type="",
    url="",
    packages=['ebspy'],
    install_requires=['matplotlib','scipy','numpy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)