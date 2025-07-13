//
//  pvdwrite.hpp
//  Version 2.0
//
//  Created by Hiroshi Sugimoto on 2018/08/28.
//  Copyright © 2018年 Hiroshi Sugimoto. All rights reserved.
//  Try It.

#ifndef pvdwrite_hpp
#define pvdwrite_hpp

#include <iostream>
#include <fstream>
#include <map>
#include <filesystem>
#include <boost/filesystem.hpp>

namespace FS=std::filesystem;

class PVDwrite {
	FS::path _filename;
	FS::path _basedir;
	using time_array=std::map<double,std::vector<FS::path>>;
	time_array _array;
	time_array::iterator find_time(double t)
	{
		for(auto itr=_array.begin();itr!=_array.end();itr++)
			if (std::abs(itr->first-t)<1e-6) return itr;
		return _array.end();
	}
public:
	PVDwrite()
	{
	};
	PVDwrite(FS::path filename):_filename(filename)
	{
		_basedir=_filename.parent_path();
	};
	~PVDwrite()
	{
		_array.clear();
	};
	PVDwrite& operator+(const PVDwrite& other)
	{
		auto part1=_filename.filename();		part1.replace_extension("");
		auto part2=other._filename.filename();	part2.replace_extension("");
		part1+="+";	part1+=part2;	part1=_basedir/part1;	part1.replace_extension("pvd");
		PVDwrite* NEW=new PVDwrite(part1);
		*NEW+=*this;
		*NEW+=other;
		return *NEW;
	}
	PVDwrite& operator+=(const PVDwrite& other)
	{
		size_t idx=0;
		for(auto& [my_time,my_name]:_array) idx=std::max(idx,my_name.size());
		for(auto& [new_time,new_name]:other._array)
		{
			auto find=find_time(new_time);
			auto& ARY=(find!=_array.end())?find->second:_array[new_time];
			ARY.resize(ARY.size()+new_name.size());
			for(size_t new_idx=0;new_idx<new_name.size();new_idx++)
				ARY[idx+new_idx]=new_name[new_idx];
		}
		return *this;
	}
	void Append(FS::path &P,double time)
	{
		_array[time].resize(1);
		_array[time][0]=P;
	};
	FS::path Filename()
	const {
		return _filename;
	}
	FS::path Filename(const FS::path& filename)
	{
		_basedir=filename.parent_path();
		return _filename=filename;
	}
	void Write()
	const {
		std::ofstream ofs(_filename);
		ofs << "<?xml version=\"1.0\"?>";
		ofs << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\" >";
		ofs << "<Collection>" << std::endl;
		for(auto& [time,name_array]:_array)
		{
			for(size_t idx=0;idx<name_array.size();idx++)
			{
				auto& name=name_array[idx];
				auto P=(_basedir.string()=="")?
					name:FS::relative(name,_basedir);
				ofs << "<DataSet timestep=\"" << time
				<< "\" group=\"\" part=\"" << idx << "\" file="
				<<	P
				<< " />" << std::endl;
			}
		}
		ofs << "</Collection></VTKFile>" << std::endl;
	};
};

#endif /* pvdwrite_hpp */
