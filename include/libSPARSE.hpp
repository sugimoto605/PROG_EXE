//
//  libSPARSE.hpp
//
//  Created by Hiroshi Sugimoto on 2018/01/14.
//  Copyright © 2018 Hiroshi Sugimoto. No rights reserved.
//

#ifndef libSPARSE_hpp
#define libSPARSE_hpp

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

namespace SPARSE
{
static double threshold=1e-14;				//これ以下は0とみなす
using index_t=std::vector<std::unordered_set<size_t>>;
class Vector  : public std::map<size_t,double>
{
public:
	void Print(size_t N,bool full=true)//N番目の添字までをプリント
	{
		int _width=10;
		for(size_t i=0;i<N;i++) {
			auto _e=find(i);
			if (_e !=end()) std::cout << "[" << i << "]" << std::fixed << std::setw(_width) << _e->second;
			else if (full) std::cout << "[" << i << "]"  << std::fixed << std::setw(_width) << "-";
		}
		std::cout<< std::endl;
	}
};
class Equation
{
	size_t _index;
	double _value;
	size_t _dimension;
	std::unordered_map<size_t,double> _less;					//小さい集合
	std::unordered_map<size_t,double> _more;					//大きい集合
public:
	Equation (size_t index,double x):_index(index),_value(x),_dimension(index+1)
	{};
	Equation(const Equation &obj)								//コピーコンストラクタ
	{
		_index=obj._index;
		_value=obj._value;
		_dimension=obj._dimension;
		for(auto &v:obj._less) _less[v.first]=v.second;
		for(auto &v:obj._more) _more[v.first]=v.second;
	}
	void AddEquation(double rate, const Equation &other, index_t &less, index_t &more)	//other式をrate倍して加える
	{
		for(auto &_j:other._less)
		{//_j.first列を自分の係数に加える
			if (_j.first<_index)
			{
				_less[_j.first]+=_j.second*rate;
				if (std::abs(_less[_j.first])>threshold*std::abs(_j.second))
					more[_j.first].insert(_index);
				else
				{
					_less.erase(_j.first);
					more[_j.first].erase(_index);
				}
			}
			else if (_j.first>_index)
			{
				_more[_j.first]+=_j.second*rate;
				if (std::abs(_more[_j.first])>threshold*std::abs(_j.second))
					less[_j.first].insert(_index);
				else
				{
					_more.erase(_j.first);
					less[_j.first].erase(_index);
				}
			}
			else
				_value+=_j.second*rate;
		}
		for(auto &_j:other._more)//_j.first列を自分の係数に加える
		{
			if (_j.first<_index)
			{
				_less[_j.first]+=_j.second*rate;
				if (_j.first<more.size())
				{
					if (std::abs(_less[_j.first])>threshold*std::abs(_j.second))
						more[_j.first].emplace(_index);
					else
					{
						_less.erase(_j.first);
						more[_j.first].erase(_index);
					}
				}
			}
			else if (_j.first>_index)
			{
				_more[_j.first]+=_j.second*rate;
				if (_j.first<less.size())
				{
					if (std::abs(_more[_j.first])>threshold*std::abs(_j.second))
						less[_j.first].emplace(_index);
					else
					{
						_more.erase(_j.first);
						less[_j.first].erase(_index);
					}
				}
			}
			else
				_value+=_j.second*rate;
		}
		{
			size_t _j=other._index;
			if (_j<_index)
			{
				_less[_j]+=other._value*rate;
				if (_j<more.size())
				{
					if (std::abs(_less[_j])>threshold*std::abs(other._value))
						more[_j].emplace(_index);
					else
					{
						_less.erase(_j);
						more[_j].erase(_index);
					}
				}
			}
			else if (_j>_index)
			{
				_more[_j]+=other._value*rate;
				if (_j<less.size())
				{
					if (std::abs(_more[_j])>threshold*std::abs(other._value))
						less[_j].emplace(_index);
					else
					{
						_more.erase(_j);
						less[_j].erase(_index);
					}
				}
			}
		}
	};
	size_t Dimension() const											//最大の添字+1
	{
		return _dimension;
	}
	const std::unordered_map<size_t,double>& map_less() const	//対角要素よりも前のインデックス
	{
		return _less;
	}
	const std::unordered_map<size_t,double>& map_more() const	//対角要素よりも後のインデックス
	{
		return _more;
	}
	bool Normalize()											//対角要素を1にする
	{
		if (std::abs(_value)>threshold)
		{
			for(auto &_e:_less) _e.second/=_value;
			for(auto &_e:_more) _e.second/=_value;
			_value=1.0;
			return true;
		}
		return false;
	};
	void Print(size_t N,bool full=true)	const					//N番目の添字までをプリント
	{
		int _width=10;
		for(size_t i=0;i<_index;i++)
		{
			auto _e=_less.find(i);
			if (_e !=_less.end())
				std::cout << "[" << i << "]" << std::fixed << std::setw(_width) << _e->second;
			else if (full)
				std::cout << "[" << i << "]"  << std::fixed << std::setw(_width) << "-";
		}
		std::cout << "[" << _index << "]" << std::fixed << std::setw(_width) << _value;
		for(size_t i=_index+1;i<N;i++)
		{
			auto _e=_more.find(i);
			if (_e !=_more.end())
				std::cout << "[" << i << "]" << std::fixed << std::setw(_width) << _e->second;
			else if (full)
				std::cout << "[" << i << "]"  << std::setw(_width) << std::fixed << "-";
		}
		std::cout<< std::endl;
	}
	const double operator[](std::size_t index) const	    	//Eq[j](方程式のj番目の係数)を引用. 存在しない係数は0を返す
	{
		if (index<_index)
		{
			auto _e=_less.find(index);
			if (_e==_less.end()) return _less.at(index);
			else return 0.;
		}
		else if (index>_index)
		{
			auto _e=_more.find(index);
			if (_e==_more.end()) return _more.at(index);
			else return 0.;
		}
		return _value;
	}
	double& operator[](std::size_t index)						//Eq[j](方程式のj番目の係数)に代入
	{
		if (index<_index)
		{
			auto _e=_less.find(index);
			if (_e==_less.end()) _less[index]=0.0;
			return _less.at(index);
		}
		else if (index>_index)
		{
			auto _e=_more.find(index);
			if (_e==_more.end()) {_more[index]=0.0; _dimension=std::max(_dimension,index+1);}
			return _more.at(index);
		}
		return _value;
	}
	void operator=(double x)                    			    //現存する係数の値を全てxにする
	{
		for(auto &_e:_less) _e.second=x;
		for(auto &_e:_more) _e.second=x;
		_value=x;
	}
	double operator*(Vector &y)							    	//係数xとベクトルyの内積を求めます.
	{
		double sum=_value*y[_index];
		for(auto &_e:_less) sum+=_e.second*y[_e.first];
		for(auto &_e:_more) sum+=_e.second*y[_e.first];
		return sum;
	}
};
class Coef
{
	std::map<size_t,Equation> _eqs;					//係数行列. 0番係数-N-1番係数からなる
	bool _index_defined=false;						//インデックスが作成されているか
	size_t _index_count;							//index
	index_t _less; 									//j列の対角要素以外では _index[j]<j行の要素が存在
	index_t _more;                                  //j列の対角要素以外では _index[j]>j行の要素が存在
	std::unordered_set<size_t> _backup;				//Erase*を行う前にはindexのバックアップが必要
	size_t _get_raw() const
	{
		size_t _raws=0;
		for(const auto &_e:_eqs) _raws=std::max(_raws,_e.second.Dimension());
		return _raws;
	}
public:
	size_t buffer_size=5;
	Coef(){};
	Coef(const Coef &obj)							//コピーコンストラクタ
	{
		_eqs=obj._eqs;
	}
	void ClearIndex()								//インデックスをクリア（削除はしない）
	{
		for(auto &_e:_less) _e.clear();
		for(auto &_e:_more) _e.clear();
		_index_defined=false;
	}
	size_t Count()									//登録されている方程式の総数を取得. Count()==Dimension()でなければ解けないよな
	{
		return _eqs.size();
	}
	size_t Dimension() const 								// 登録されている方程式のIDの最大値を取得. 要するに係数行列の行数
	{
		//方程式系の中で, 最大の対角要素のindex+1が方程式の数である.
		if (_eqs.size()==0) return 0;
		size_t ret=0;
		for(auto &_e:_eqs) {
			ret=std::max(ret,_e.first);
		}
		return ret+1;								// 方程式総数が2で，ID=2, 5なら, Dimension()=6である. ID=0,1,3,4の方程式の追加必要
	};
	bool EraseDown()
	{
		for(size_t j=0;j<_index_count-1;j++) if (!EraseDown(j)) return false;
		return true;
	}
	bool Swap(size_t i,size_t j)					// i行目とj行目を交換する
	{
		if (!_index_defined) throw(std::runtime_error("libSPARSE::Swap::Can not find index."));
		if (i==j) return true;						//同じ行ならば何もしない
		auto i_eq=_eqs.extract(i);					//i行目の方程式を取り出す
		auto j_eq=_eqs.extract(j);					//j行目の方程式を取り出す
		i_eq.key()=j;								//i行目の方程式のIDをjに変更
		j_eq.key()=i;								//j行目の方程式のIDをiに変更
		_eqs.insert(std::move(i_eq));				//i行目の方程式を挿入
		_eqs.insert(std::move(j_eq));				//j行目の方程式を挿入
		for(auto &m:_more)
		{
			m.erase(i);
			m.insert(j);
			m.erase(j);
			m.insert(i);
		}
		for(auto &m:_less)
		{
			m.erase(i);
			m.insert(j);
			m.erase(j);
			m.insert(i);
		}
		return true;
	};
	bool EraseDown(size_t j)						//j行を用いて, j列の下三角部分を消去
	{
		if (!_index_defined) throw(std::runtime_error("indexが未完成です."));
		//j,j+1,...,_index_count-1行を操作する. j列目の絶対値が一番大きいものを, j列目に交換する
		size_t pivot=j;									//j行目をピボットとする
		double max_value=std::abs(_eqs.at(pivot)[j]);
		for(auto i:_more[j]) if (std::abs(_eqs.at(i)[j])>max_value) max_value=std::abs(_eqs.at(pivot=i)[j]);
		if (std::abs(max_value) < threshold) return false;	//ピボットが0ならば消去できない
		// if (j != pivot)
		// {
		// 	std::cout << "replace line " << j << " with should be replaced by pivot=" << pivot << std::endl;
		// 	auto j_eq=_eqs.extract(j);		//j行目の方程式を取り出す
		// 	auto p_eq=_eqs.extract(pivot);	//pivot行目の方程式を取り出す
		// 	j_eq.key()=pivot;				//j行目の方程式のIDをpivotに変更
		// 	p_eq.key()=j;					//pivot行目の方程式のIDをjに変更
		// 	_eqs.insert(std::move(j_eq));	//j行目の方程式を挿入
		// 	_eqs.insert(std::move(p_eq));	//pivot行目の方程式を挿入
		// 	//インデックスを更新
		// 	for(auto& m:_more)
		// 	{

		// 	}
		// 	for (size_t k = 0; j < _index_count - 1; j++)
		// 		for (auto i : _more[j])
		// 		{
		// 			_less[i].erase(j);
		// 			_less[i].insert(pivot);
		// 		}
		// }
		//AddRquation()がインデックスを書き換えるのでコピーしておく
		_backup.clear();
		for(auto i:_more[j]) _backup.insert(i);
		//j行を用いてi行を削除 (_eqs[j]だとconstにならずコピーが生じる)
		for(auto i:_backup)
		{
			if (!_eqs.at(j).Normalize()) return false;
			_eqs.at(i).AddEquation(-_eqs.at(i)[j]/_eqs.at(j)[j], _eqs.at(j),_less,_more);
		}
		_eqs.at(_index_count-1).Normalize();
		return true;
	};
	bool EraseUp()
	{
		for(size_t j=_index_count-1;j>0;j--) if (!EraseUp(j)) return false;
		return true;
	}
	bool EraseUp(size_t j)							//j行を用いて, j列の上三角部分を消去
	{
		if (!_index_defined) throw(std::runtime_error("indexが未完成です."));
		//インデックスを書き換える前にコピー
		_backup.clear();
		for(auto i:_less[j]) _backup.insert(i);
		//j行を用いてi行を削除 (_eqs[j]だとconstにならずコピーが生じる)
		for(auto i:_backup)
		{
			if (!_eqs.at(j).Normalize()) return false;
			_eqs.at(i).AddEquation(-_eqs.at(i)[j]/_eqs.at(j)[j], _eqs.at(j),_less,_more);
		}
		_eqs.at(0).Normalize();
		return true;
	};
	void RemoveIndex()								//インデックスを削除(再作成)
	{
		while (!_less.empty())
		{
			_less.back().clear();
			_less.pop_back();
		}
		while (!_more.empty())
		{
			_more.back().clear();
			_more.pop_back();
		}
		_index_defined=false;
	}
	bool CheckIndex() const						//インデックスが正しいかチェック
	{
		if (!_index_defined) 
		{
			std::cerr << "libSPARCE::CheckIndex::Can not find index." << std::endl;
			return false;
		}
		for(size_t i=0;i<_index_count;i++)
		{
			for (auto &_e : _less[i])
			{
				std::cout << "CheckIndex.less[" << i << "][" << _e << "]" << std::endl;
				if (_e >= i)
				{
					std::cerr << "less["<<i<<"]element[" << _e << "]must be smaller than[" << i << "]" << std::endl;
					return false; // less集合はi行目より小さい
				}
				if (std::abs(_eqs.at(_e)[i])<threshold)
				{
					std::cerr << "Matrix[" << _e << "][" << i << "]= " << _eqs.at(_e)[i] << " close to zero." << std::endl;
				}
			}
			for(auto &_e:_more[i])
			{
				std::cout << "CheckIndex.more[" << i << "][" << _e << "]" << std::endl;
				if (_e <= i)
				{
					std::cerr << "more["<<i<<"]element[" << _e << "]must be larger than[" << i << "]" << std::endl;
					return false; // more集合はi行目より大きい
				}
				_eqs.at(_e).Print(_get_raw());
				for(size_t kk=0;kk<_get_raw();kk++)
				{
					std::cout << "VAL[" << kk << "]=" << _eqs.at(_e)[kk] << std::endl;
				}
				if (std::abs(_eqs.at(_e)[i])<threshold)
				{
					std::cerr << "Matrix[" << _e << "][" << i << "]= " << _eqs.at(_e)[i] << " is close to zero." << std::endl;
				}
			}
		}
		return true;
	}
	void MakeIndex()								//インデックスを作成
	{
		if (_index_defined) ClearIndex();
		_index_count=Dimension();
		_less.reserve(_index_count);
		for(size_t i=_less.size();i<_index_count;i++) _less.emplace_back();
		for(auto &_e:_less) _e.reserve(buffer_size);               //普通これだけあれば大丈夫
		_more.reserve(_index_count);
		for(size_t i=_more.size();i<_index_count;i++) _more.emplace_back();
		for(auto &_e:_more) _e.reserve(buffer_size);               //普通これだけあれば大丈夫
		for(auto &_e:_eqs)
		{
			for(auto &_f:_e.second.map_more()) //上三角部：i行=_e.first, j列=_f.firstの値が_f.second
				if (std::abs(_f.second)>threshold)
					if (_f.first<_index_count)
						_less[_f.first].insert(_e.first);
			for(auto &_f:_e.second.map_less()) //下三角部：i行=_e.first, j列=_f.firstの値が_f.second
				if (std::abs(_f.second)>threshold)
					if (_f.first<_index_count)
						_more[_f.first].insert(_e.first);
		}
		_backup.reserve(buffer_size);
		_index_defined=true;
	}
	void New(size_t index, double x)				//対角要素xの方程式を追加
	{
		if (_index_defined) throw(std::runtime_error("index作成後は不可です"));
		auto _e=_eqs.find(index);
		if (_e!=_eqs.end()) throw(std::runtime_error("その方程式, もうあるやんけ"));
		_eqs.emplace(std::piecewise_construct,std::forward_as_tuple(index),std::forward_as_tuple(index,x));
	}
	void Print(bool full=true)						//係数行列の次元まで全てプリント
	{
		size_t lines=Dimension();//	Dimension()は行数
		size_t raws=_get_raw();
		std::cout << lines << "x" << raws << std::endl;
		for(size_t i=0;i<lines;i++) {
			std::cout << "(" << std::setw(4)  << std::setfill(' ') << i << ") ";
			auto _eq=_eqs.find(i);
			if (_eq!=_eqs.end()) {
				_eq->second.Print(raws,true);
			} else if (full) std::cout <<  std::endl ;
		}
	}
	void PrintIndex()								//インデックスを表示
	{
		//indexは, 行列の各行について定義されている.
		//i行目のindexとは, i列目が非ゼロである行の集合である.
		//		0..i-1行目の集合をless集合 _less[i] と呼び, i+1...N-1行目の集合をmore集合 _more[i] と呼ぶ.
		//	    集合は, std::unordered_set<size_t>である.
		if (_index_defined)
		{
			std::cout << "idx=";
			for(int j=0; j<_index_count;j++) {
				std::cout << "[" << j << ",less";
				for(auto &_e:_less[j]) std::cout << ":" << _e ;
				std::cout << ",more";
				for(auto &_e:_more[j]) std::cout << ":" << _e ;
				std::cout << "]";
			}
			std::cout << std::endl;
		}
	}
	const Equation& operator[](std::size_t index) const //i行目の係数を引用
	{
		return _eqs.at(index);
	};
	Equation& operator[](std::size_t index)			//i行目の係数を設定
	{
		if (_index_defined) throw(std::runtime_error("index作成後は不可です"));
		auto _e=_eqs.find(index);
		if (_e==_eqs.end()) _eqs.emplace(std::piecewise_construct,std::forward_as_tuple(index),std::forward_as_tuple(index,0.0));
		return _eqs.at(index);
	};
	void operator=(double x)						//現存する係数の値をxにする
	{
		if (_index_defined) throw(std::runtime_error("index作成後は不可です"));
		for(auto &_e:_eqs) _e.second=x;
	};
	Vector& operator*(Vector &x)					//係数行列にxを掛ける
	{
		Vector *_y=new Vector();
		for(auto &_e:_eqs) (*_y)[_e.first]=_e.second*x;
		return *_y;
	};
};

}//end namespace SPARSE

//void test_libSPARSE()
//{
//	std::cout << "--Vectorクラスで疎ベクトルを保存--" << std::endl;
//	SPARSE::Vector vec;
//	vec[0]=2.0;	vec[4]=1.0;
//	vec.Print(6);
//	std::cout << "--Equationクラスで疎係数を保存--" << std::endl;
//	SPARSE::Equation eq(3,3.0);
//	eq[5]=3.2;
//	eq.Print(6);
//	std::cout << "  Equation.Dimension=" << eq.Dimension() << std::endl;
//	std::cout << "  0でクリアできるよ--" << std::endl;
//	eq=0;
//	eq.Print(6);
//	std::cout << "  ベクトルや係数にスカラー実数を代入すると, かつて成分が存在したところだけに代入される謎仕様" << std::endl;
//	eq=2;
//	eq.Print(6);
//	//-----------------------
//	std::cout << std::endl << std::endl
//	<< "まあおふざけはこの程度にしておいて，疎な行列のソルバー作る" << std::endl;
//	std::cout << "プログラミング効率を優先し, 実行性能を犠牲にしたものだ" << std::endl;
//	std::cout << "--Coefクラスで係数行列を保存--" << std::endl;
//	SPARSE::Coef coef;
//	for(int i=0;i<6;i++) coef.New(i,1.0); // 6x6 のときは6回New()する
//	coef[0][2]=-1; coef[0][1]=2.3; coef[0][2]=-1;
//	coef[0][3]=+1; coef[0][4]=-3; coef[0][5]=+1;
//	coef[1][0]=-2.3; coef[1][2]=3.0;
//	coef[2][1]=-3.0; coef[2][3]=3.0;
//	coef[3][2]=-3.0; coef[3][4]=3.0;
//	coef[4][3]=-3.0; coef[4][5]=3.0;
//	coef[5][4]=-3.0; coef[5][0]=-1.0;
//	coef.Print();   //変数が存在しないところは - と表示される
//	std::cout << "--coef x Vector は 行列 x ベクトル--" << std::endl;
//	auto res=coef*vec;
//	res.Print(6);
//	std::cout << "  なお, 積を行うと疎ベクトルが密ベクトルに変化することに注意:" << std::endl;
//	vec.Print(6);
//	std::cout << std::endl << "--係数Coef の末尾に非斉次項を追加して連立方程式を定義する--"  << std::endl;
//	for(int i=0;i<coef.Dimension();i++) coef[i][coef.Dimension()]=res[i];
//	std::cout << "  なお, 係数Coef の末尾に, 好きにおまけ領域を確保しても良い. 旅も計算も道連れだ"  << std::endl;
//	for(int i=0;i<coef.Dimension();i++) coef[i][coef.Dimension()+1]=i*i;
//	coef.Print();
//	std::cout << "--最初にindexを作成" << std::endl;
//	coef.MakeIndex();
//	coef.PrintIndex();
//	std::cout << "--下三角行列に変形" << std::endl;
//	coef.EraseUp();
//	coef.Print();
//	coef.PrintIndex();
//	std::cout << "--上三角行列に変形" << std::endl;
//	coef.EraseDown();
//	coef.Print();
//	coef.PrintIndex();
//	std::cout << "--ロック解除--" << std::endl;
//	coef.ClearIndex();
//	std::cout << "--解は非斉次部分に並ぶ--" << std::endl;
//	SPARSE::Vector solution;
//	for(int i=0;i<coef.Dimension();i++) solution[i]=coef[i][coef.Dimension()];
//	solution.Print(coef.Dimension());
//	std::cout << "--念のため係数消しとこ--" << std::endl;
//	coef=0;
//	coef.Print();
//	coef.PrintIndex();
//}


#endif /* libSPARSE_hpp */
