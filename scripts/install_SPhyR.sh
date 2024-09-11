SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# install LEMON first
cd $SCRIPT_DIR/../oncophylo/lib
echo "Downloading LEMON in oncophylo/lib"
curl "http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz" --output lemon-1.3.1.tar.gz

echo "Installing LEMON in oncophylo/lib"
tar xvzf lemon-1.3.1.tar.gz
cd lemon-1.3.1
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$SCRIPT_DIR/../oncophylo/lib/lemon
make
make check
make install

# install SPhyR
echo "Installing SPhyR in oncophylo/lib"
cd $SCRIPT_DIR/../oncophylo/lib/SPhyR
mkdir build
cd build
cmake .. -DLIBLEMON_ROOT=$SCRIPT_DIR/../oncophylo/lib/lemon
make

mv $SCRIPT_DIR/../oncophylo/lib/SPhyR/build/kDPFC $SCRIPT_DIR/../oncophylo/bin
mv $SCRIPT_DIR/../oncophylo/lib/SPhyR/build/visualize $SCRIPT_DIR/../oncophylo/bin

echo "IF THIS INSTALL FAILS BECAUSE OF THE FOLLOWING ERROR: "
echo "          ~/OncoPhylo/oncophylo/lib/SPhyR/src/columngen.h:11:10: fatal error: 'ilcplex/ilocplex.h' file not found 
                #include <ilcplex/ilocplex.h>"
echo "You'll need to make sure CPLEX is install and that the CPLEX_AUTODETECTION_PREFIXES in CMakeLists.txt contains the directory where CPLEX is install on your machine."