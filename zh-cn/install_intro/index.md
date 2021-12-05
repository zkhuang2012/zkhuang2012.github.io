# 关于生信软件的安装


生信学习里装软件是一件很头疼又富有成就感的事，且整理一下一些心得。

## 环境部署

开始之前，建议一开始先把工作目录都规划好(养成好习惯)。

```{shell}
# 部署工作目录
mkdir -p {bin,download,local/{app,include,lib},project,perl5}

tree --charset ASCII -L 2

# 如下结构
.
|-- bin
|-- download
|-- local
|   |-- app
|   |-- include
|   `-- lib
|-- perl5
|-- project

##bin 直接添加到PATH路径中，后续安装完软件，直接```ln -s <path-to-soft> $HOME/bin/```
##down  用来下载软件和数据常用,下载命令如**wget,axel**等
##local 用来装软件和各种库，app安装软件
##pel5  用来装perl module，后面会介绍
##project 用来管理和分析项目
```

建好文件夹后就要在开始进行环境变量设置，即在**$HOME/.bashrc**中添加如下内容，添加后记得```source $HOME/.bashrc```

```
# Local bin,app
export LOCAL_BIN=$HOME/bin
export LOCAL_HOME=$HOME/local
export LOCAL_APP=$LOCAL_HOME/app
export PATH=$LOCAL_BIN:$PATH

# Local lib,include
## 动态链接库的路径
export LD_LIBRARY_PATH=$LOCAL_HOME/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LOCAL_HOME/lib:$LD_RUN_PATH
## 静态链接库的路径
export LIBRARY_PATH=$LOCAL_HOME/lib:$LIBRARY_PATH
## gcc头文件的路径
export C_INCLUDE_PATH=$LOCAL_HOME/include
## g++头文件的路径
export CPLUS_INCLUDE_PATH=$LOCAL_HOME/include:$CPLUS_INCLUDE_PATH

# Local perl lib
LOCAL_PERL_EDITION=${HOME}/perl5
export PERL5LIB=$LOCAL_PERL_EDITION/lib
export PATH=$LOCAL_PERL_EDITION/bin:$PATH

```
以上部署时N年前参考[shenwei博客](http://blog.shenwei.me/about/)等，感兴趣可以去瞅瞅。

## 软件安装
**1st**: 最希望的是软件开发者都编译好了，这样我们直接下载解压后就能用那种的，比如blast+等。

**2nd**: 源代码(source code)安装，解压后，一定要好好看README.txt或者INSTALL，一般来说就是三部曲。
```
# 1.检测环境依赖库
./configure --prefix=/path-to-install/
## configure中通过--prefix=/path-to-install/来指定安装路径

# 2.编译软件
make -j4

# 3.安装软件
make install
```
- Q1:configure过程往往都是各种error或者not found等错误
- A1:这样需要通过错误提示去搜索，可能缺某个lib，一般解决方法就是下载lib source code安装，然后把lib加入环境的动态或者静态库路径中。
- Q2:有的程序没有configure直接make，那怎么指定安装路径呢？
- A2:这样可以通过make install DESTDIR=/path-to-install/来指定安装路径


**3rd**: 这时候还有conda，因为系统默认是python2，所以下个python2版本的conda
```
# 创建环境，如env_test
conda create -n env_test

# 激活环境
conda activate env_test

# 安装软件
conda install python=3.6
##conda install <软件名>=<version>
```

至于源代码安装和conda安装的优先级，我个人倾向源代码，虽然要经常debug。
conda有以下两点，
- 依赖镜像地址，一旦地址失效，就只能离线安装；
- 环境设置不好容易混乱。

总之，仁者见仁，能用就行。

接下来，以实验室集群的环境进行案例演示

### 案例一：R语言安装
如果你的工作目录部署和环境配置是和我上面介绍的一样的话，那么，这下面R-3.6.3的安装代码**直接复制到我们的集群上运行**即可。
```
# 建立软件目录
mkdir -p $LOCAL_APP/R
## 经常要装多个R版本，所以索性建个R文件夹，子文件夹再分版本

# 下载软件
wget -c $HOME/download https://mirrors.e-ducation.cn/CRAN/src/base/R-3/R-3.6.3.tar.gz

# 解压软件
tar -xvzf $HOME/download/R-3.6.3.tar.gz -C $LOCAL_APP/R

# 指定GCC环境，集群公共软件里已经有安装gcc-7.2,所以直接设置添加到.bashrc中

cat <<EOF >>$HOME/.bashrc

# gcc
LOCAL_GCC_PATH=/public/tools/devtools/gcc-7.2.0/
export PATH=\$LOCAL_GCC_PATH/bin:\$PATH
export C_INCLUDE_PATH=\$LOCAL_GCC_PATH/lib:\$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=\$LOCAL_GCC_PATH/lib64:\$CPLUS_INCLUDE_PATH
export CFLAGS=-I\$LOCAL_GCC_PATH/include:\$CFLAGS
export LDFLAGS=-L\$LOCAL_GCC_PATH/lib64:-L\$LOCAL_GCC_PATH/lib:\$LDFLAGS
export LD_LIBRARY_PATH=\$LOCAL_GCC_PATH/lib64:\$LD_LIBRARY_PATH
export LD_RUN_PATH=\$LOCAL_GCC_PATH/lib64:\$LD_RUN_PATH
export LIBRARY_PATH=\$LOCAL_GCC_PATH/lib64:\$LIBRARY_PATH

EOF

source $HOME/.bashrc

# 开始安装，这里我默认安装在R语言安装目录，所以用--prefix=$PWD,大家也可以修改--prefix=/path-to-install-R/

cd $LOCAL_APP/R/R-3.6.3

./configure --prefix=$PWD --enable-R-shlib --with-blas --with-lapack --with-x  --with-tcltk --with-tcl-config=/usr/lib64/tcl8.5/tclConfig.sh --with-tk-config=/usr/lib64/tkConfig.sh

make -j4 && make install

# 添加到路径
ln -s $PWD/bin/R $LOCAL_BIN/R
ln -s $PWD/bin/Rscript $LOCAL_BIN/Rscript

# 配置R语言环境,即修改$HOME/.Rprofile文件
cat <<EOF >$HOME/.Rprofile

local({
  r = getOption("repos")
  r["CRAN"] = "http://mirrors.ustc.edu.cn/CRAN/"
  options(repos = r)
})
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

EOF

# END
```

### 案例二：Python安装
-	Q1:首先，要装python2还是python3？
-	A1:一般是后者，因为集群已经自带python 2.7.15，且python3是大趋势。
-	Q2:其次，虽然系统自带python2，但普通用户没有权限将module安装到系统Python库，怎么破？
-	A2:两种方法，一种就是通过安装命令参数和环境配置解决，一种源代码安装（案例装的python3）。

#### 通过命令参数和环境变量设置安装私有包
```
# 安装私有包
##通过pip安装，简单粗暴
##pip install --user <packagename>
##下载module source code然后执行以下命令
##python3 setup.py install --user

# 举例装个matplotlib，这里-i加入清华镜像，可以提速
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple matplotlib


# 安装好的私有包就在$HOME/.local/lib/python2.7/site-packages,添加到环境变量即可
cat <<EOF >>$HOME/.bashrc

# python2.7 module path
export PYTHONPATH=\$HOME/.local/lib/python2.7/site-packages:\$PYTHONPATH

EOF

source $HOME/.bashrc
```

#### 源代码安装python3
如果系统连python2安装缺少某些库，那一样只能源代码安装python，这里以安装python3.6为案例。

```
# 下载软件安装包
axel -a -o $HOME/download https://www.python.org/ftp/python/3.6.10/Python-3.6.10.tgz
# 解压压缩包
tar -xvzf $HOME/download/Python-3.6.10.tgz -C $LOCAL_APP/
# 安装
./configure --prefix=$PWD --enable-shared --with-ensurepip=install --with-tcltk-includes="-I/usr/include/" --with-tcltk-libs="-L/usr/lib64/"
make -j4 && make install

cat <<EOF >>$HOME/.bashrc

# local python3 path
export LOCAL_PYTHON=\$LOCAL_APP/Python-3.6.10
export PYTHONPATH=\$LOCAL_PYTHON/lib/python3.6/site-packages:$LOCAL_PYTHON/lib:$PYTHONPATH
export PATH=\$LOCAL_PYTHON/bin:\$PATH
export LD_LIBRARY_PATH=\$PYTHONPATH/lib:\$LD_LIBRARY_PATH
export LD_RUN_PATH=\$LOCAL_GCC_PATH/lib64:\$LD_RUN_PATH

EOF

source $HOME/.bashrc
```
PS：python2和python3是相互不兼容的额，所以当装了多个python，在切换版本时，一定要确认$PYTHONPATH是否也对应切换了。

### 案例三：Perl module安装
既然python可以利用系统python，perl当然也可以。
但不建议自己安装perl，过程比较麻烦，还要注意安装perl版本与系统版本是一致的，不然有的软件安装默认调用系统perl，一旦模块对不同版本不兼容，就会各种报错，相当郁闷。
```
# 下载cpanm
curl -L http://cpanmin.us/ -o $LOCAL_BIN/cpanm && chmod +x $LOCAL_BIN/cpanm
##cpanm是安装Perl模块的最方便的方法，自动下载安装依赖包。
##cpanm其实是一个可执行文件而已，将它下载到bin目录，然后添加执行权限就可以了。

# 设置环境变量
cat <<EOF >>$HOME/.bashrc

# local perl lib path
alias  cpanm='cpanm --prompt --mirror http://mirrors.ustc.edu.cn/CPAN/ -l  \$HOME/perl5'
export PERL5LIB=\$HOME/perl5/lib/perl5:\$PERL5LIB

EOF

source $HOME/.bashrc

# 安装模块命令cpanm <module-name>，这里举例安装MongoDB
cpanm MongoDB

# 如果要卸载模块则要借助App::pmuninstal（也是一个模块）
cpanm App::pmuninstal
# 卸载模块命令pm-uninstal <module-name>
pm-uninstal MongoDB

```

## 参考
-	http://blog.shenwei.me/about/
-	http://www.ttlsa.com/perl/use-cpanm-to-install-perl-modules/

## END

