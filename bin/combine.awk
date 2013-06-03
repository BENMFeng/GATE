BEGIN{
	i=0;
	c="";
	OFS="\t";
}
{
	if (c!=$1){
		if (c != ""){
			combine(chr,c,i);
		}
		i=0;
		c="";
	}
	c=$1;
	chr[i,0]=$2;
	chr[i,1]=$3;
	i++
}
END{
	combine(chr,c,i);
}
function combine(arr,c,i) {
	j=0;
	new[j,0]=arr[0,0];
	new[j,1]=arr[0,1];
	for (k=1;k<i;k++)
	{
		if ((arr[k,0]<=new[j,1])&&(arr[k,1]>=new[j,1])){
			new[j,1]=arr[k,1];
		}
		else if (arr[k,0]>new[j,1]){
			j++;
			new[j,0]=arr[k,0];
			new[j,1]=arr[k,1];
		}
	}
	for (n=0;n<=j;n++){
		print c,new[n,0],new[n,1]
	}
}