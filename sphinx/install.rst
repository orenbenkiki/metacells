Installation
============

In short: ``pip install metacells``. If you do not have ``sudo`` privileges, you might need to ``pip
install --user metacells``. Note that ``metacells`` requires many "heavy" dependencies, most notably
``numpy``, ``pandas``, ``scipy``, ``scanpy``, ``leidenalg``, which ``pip`` should automatically
install for you.

The metacells package contains extensions written in C++, which means installing it requires
compilation. In Ubuntu the necessary tools can be installed using ``sudo apt-get install
build-essential`` - please refer to the documentation for your specific system. Make sure that not
only ``gcc`` is installed, but also ``g++`` to allow for C++ compilation.

.. todo::

    Provide pre-build wheel binaries, especially for the unfortunates using Windows.

You can check if all went well by running a variety of tests using ``pytest`` or ``tox``.
