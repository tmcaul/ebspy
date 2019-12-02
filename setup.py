import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='ebspy',
    version="1.0.0",
    author="Tom McAuliffe",
    author_email="tmcauliffe@icloud.com",
    description="A package for handling electronbackscatter diffraction data",
    long_description=long_description,
    long_description_content_type="",
    url="",
    packages=['ebspy'],
    install_requires=['matplotlib','scipy','numpy','pillow','h5py','pathlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache 2.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)