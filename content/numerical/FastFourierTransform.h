/**
 * Author: Teetat T.
 * Date: 2024-01-19
 * Description: Fast Fourier transform
 * Time: $O(N \log N)$
 */
#pragma once

namespace fft{

using cd = complex<db>;
using vc = vector<cd> ;

void fft(vc &a){
  	int n=a.size(),L=31-__builtin_clz(n);
  	vc rt(n);
  	rt[1]=1;
  	for(int k=2;k<n;k*=2){
    	cd z=polar(db(1),PI/k);
    	for(int i=k;i<2*k;i++)rt[i]=i&1?rt[i/2]*z:rt[i/2];
  	}
	vi rev(n);
	for(int i=1;i<n;i++)rev[i]=(rev[i/2]|(i&1)<<L)/2;
	for(int i=1;i<n;i++)if(i<rev[i])swap(a[i],a[rev[i]]);
	for(int k=1;k<n;k*=2)for(int i=0;i<n;i+=2*k)for(int j=0;j<k;j++){
		cd z=rt[j+k]*a[i+j+k];
		a[i+j+k]=a[i+j]-z;
		a[i+j]+=z;
	}
}

template<class T>
struct _conv{
	using vt = vector<T>;
	static const bool INT=true;
	template<class U>
	db norm(const U &x){
		return INT?round(x):x;
	}
	vt operator()(const vt &a,const vt &b){
		if(a.empty()||b.empty())return {};
		vt res(a.size()+b.size()-1);
		int L=32-__builtin_clz(res.size()),n=1<<L;
		vc in(n),out(n);
		copy(a.begin(),a.end(),in.begin());
		for(int i=0;i<b.size();i++)in[i].imag(b[i]);
		fft(in);
		for(auto &x:in)x*=x;
		for(int i=0;i<n;i++)out[i]=in[-i&(n-1)]-conj(in[i]);
		fft(out);
		for(int i=0;i<res.size();i++)res[i]=norm(imag(out[i])/(4*n));
		return res;
	}
};
template<>
struct _conv<db>{
	static const bool INT=false;
};
template<class T>
vector<T> conv(const vector<T> &a,const vector<T> &b){
	return _conv<T>()(a,b);
}

template<ll M=MOD>
struct _convMod{
	vl operator()(const vl &a,const vl &b){
		if(a.empty()||b.empty())return {};
		vl res(a.size()+b.size()-1);
		int L=32-__builtin_clz(res.size()),n=1<<L;
		ll cut=int(sqrt(M));
		vc in1(n),in2(n),out1(n),out2(n);
		for(int i=0;i<a.size();i++)in1[i]=cd(ll(a[i])/cut,ll(a[i])%cut); // a1 + i * a2
		for(int i=0;i<b.size();i++)in2[i]=cd(ll(b[i])/cut,ll(b[i])%cut); // b1 + i * b2
		fft(in1),fft(in2);
		for(int i=0;i<n;i++){
			int j=-i&(n-1);
			out1[j]=(in1[i]+conj(in1[j]))*in2[i]/(2.l*n); // f1 * (g1 + i * g2) = f1 * g1 + i f1 * g2
			out2[j]=(in1[i]-conj(in1[j]))*in2[i]/cd(0.l,2.l*n); // f2 * (g1 + i * g2) = f2 * g1 + i f2 * g2
		}
		fft(out1),fft(out2);
		for(int i=0;i<res.size();i++){
			ll x=round(real(out1[i])),y=round(imag(out1[i]))+round(real(out2[i])),z=round(imag(out2[i]));
			res[i]=((x%M*cut+y)%M*cut+z)%M; // a1 * b1 * cut^2 + (a1 * b2 + a2 * b1) * cut + a2 * b2
		}
		return res;
	}
};
template<int M=MOD>
vl conv(const vl &a,const vl &b){
	return _convMod<M>()(a,b);
}

} // namespace fft
