#!/bin/bash
#Record the beginning time
date
declare -i count=0

#Set the two-way coupled run dir and executable
#rundir="template_dir_to_be_replaced_automatically_later"
rundir='/data/users/yanyy/v11/twoway_test/'
dir_nested_NA="NA/"
dir_nested_CH="CH/"
dir_nested_EU="EU/"
dir_2x25="2x25/"
dir_4x5="4x5/"


mkdir $rundir/lock/
rm -f $rundir/lock/*

mkdir $rundir/exchange/
rm -f $rundir/exchange/*

#Options for two-way coupled runs
while getopts ":24cne" optname
do
	case $optname in
		c)
			if [ "$count" == 0 ];then
				nested1="NESTED_CH"
				nested1_dir=$rundir$dir_nested_CH
				count=1
			elif [ "$count" == 1 ];then
				nested2="NESTED_CH"
				nested2_dir=$rundir$dir_nested_CH
				count=2
			elif [ "$count" == 2 ];then
				nested3="NESTED_CH"
				nested3_dir=$rundir$dir_nested_CH
				count=3
			fi
			;;
		n)
			if [ "$count" == 0 ];then
				nested1="NESTED_NA"
				nested1_dir=$rundir$dir_nested_NA
				count=1
			elif [ "$count" == 1 ];then
				nested2="NESTED_NA"
				nested2_dir=$rundir$dir_nested_NA
				count=2
			elif [ "$count" == 2 ];then
				nested3="NESTED_NA"
				nested3_dir=$rundir$dir_nested_NA
				count=3
			fi
			;;
                e)
			if [ "$count" == 0 ];then
				nested1="NESTED_EU"
                                nested1_dir=$rundir$dir_nested_EU
				count=1
			elif [ "$count" == 1 ];then
				nested2="NESTED_EU"
				nested2_dir=$rundir$dir_nested_EU
				count=2
			elif [ "$count" == 2 ];then
				nested3="NESTED_EU"
				nested3_dir=$rundir$dir_nested_EU
				count=3
			fi
			;;
		2)
			global_dir=$rundir$dir_2x25
			;;
		4)
			global_dir=$rundir$dir_4x5
			;;
		*)
			echo "unknown arguement"
			echo 'usage : ./run.sh [ -c ] [ -n ] [ -e ] '
			exit 1
			;;
	esac
done

echo "                 GLOBAL is exchanging with $nested1 $nested2 $nested3"

pushd $global_dir
echo "  @dir: $global_dir"
./geos > log.geos &
popd

if [ "$count" == 1 ] || [ "$count" == 2 ] || [ "$count" == 3 ];then
	pushd $nested1_dir
	echo "  @dir: $nested1_dir"
	./geos > log.geos &
	popd
fi
if [ "$count" == 2 ] || [ "$count" == 3 ];then
	pushd $nested2_dir
	echo "  @dir:$nested2_dir"
	./geos > log.geos &
	popd
fi
if [ "$count" == 3 ];then
	pushd $nested3_dir
	echo "  @dir: $nested3_dir"
	./geos > log.geos &
	popd
fi

