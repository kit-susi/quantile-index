wget ftp://ftp.kernel.org/pub/linux/kernel/v2.6/linux-2.6.11.6.tar.gz
wget ftp://ftp.gnu.org/gnu/gcc/gcc-4.0.0/gcc-4.0.0.tar.bz2

tar xvzf linux-2.6.11.6.tar.gz
tar xvfj gcc-4.0.0.tar.bz2

touch sources
find linux-2.6.11.6 -regex '.*\.\(c\|h\|C\|java\)$' | while read line; do
    cat $line >> sources
    echo "\1" >> sources 
done
find gcc-4.0.0 -regex '.*\.\(c\|h\|C\|java\)$' | while read line; do
    cat $line >> sources
    echo "\1" >> sources 
done
./../../build/release/convert2surf sources
rm linux-2.6.11.6.tar.gz gcc-4.0.0.tar.bz2 
rm linux-2.6.11.6 gcc-4.0.0 -r
rm sources
echo "Done"
