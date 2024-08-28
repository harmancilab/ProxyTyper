# Installing ProxyTyper

ProxyTyper contains an executable (*bin/ProxyTyper_Release*) and a bash script (*scripts/ProxyTyper.sh*) that implements the user-facing options using the executable.

## Clone the Repository
First, we need to clone the git repository:
```
git clone https://github.com/harmancilab/ProxyTyper.git
```



## Building *ProxyTyper_Release*
To build the executable, run the following commands under ProxyTyper's main directory:
```
cd ProxyTyper
make clean
make
```
this command builds the executable and places it under */bin/ProxyTyper_Release*. Building the executable requires *zlib* on the system.

In order to run *ProxyTyper.sh* from any directory, it is necessary to add *ProxyTyper_Release* to your path.

Following command modifies .bash_profile to add executable directory to the path.
```
echo "export PATH=\$PATH:$PWD/bin" >> $HOME/.bash_profile
source $HOME/.bash_profile
```

After running these, test if *ProxyTyper_Release* is in your path:
```
ProxyTyper_Release
```
this command should write numerous options and exit. 

If you have root privileges, you can also copy *ProxyTyper_Release* to a directory in your path, e.g., */usr/bin*.

## Running *ProxyTyper.sh*
After executable is built and installed, you can copy *ProxyTyper.sh* under the directory that contains the VCF file(s). *ProxyTyper.sh* requires the configuration file named *PROXYTYPER.ini* to be present in the same directory. After copying *ProxyTyper.sh* and *PROXYTYPER.ini*, test the script by:
```
./ProxyTyper.sh
```
which should show options to run *ProxyTyper.sh* and exit with no errors. 

---
After installation, we recommend reviewing examples under *\examples* folder.