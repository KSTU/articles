#using Plots
using LsqFit
using Gadfly
using Cairo
using LaTeXStrings
#using LaTeXString

#gr()
file_id=open("sc","r")
lines=0
for i in eachline(file_id)
	global lines+=1
end
close(file_id)
print("$(lines)")

sc_dens=Array{Float64}(undef,lines)
sc_xa=Array{Float64}(undef,lines)
sc_da=Array{Float64}(undef,lines)
sc_dda=Array{Float64}(undef,lines)
sc_db=Array{Float64}(undef,lines)
sc_ddb=Array{Float64}(undef,lines)

file_id=open("sc","r")
id=0
for i in eachline(file_id)
#	tempstring=readline(file_id)
	global id+=1
#	temp=split(i)[1]
	sc_dens[id]=parse(Float64,split(i)[1])
#	println("$(temp) ||")
	sc_xa[id]=parse(Float64,split(i)[2])
	sc_da[id]=parse(Float64,split(i)[3])
	sc_dda[id]=parse(Float64,split(i)[4])
	sc_db[id]=parse(Float64,split(i)[5])
	sc_ddb[id]=parse(Float64,split(i)[6])
end
close(file_id)

######################33
file_id=open("bcc","r")
lines=0
for i in eachline(file_id)
	global lines+=1
end
close(file_id)
print("$(lines)")

bcc_dens=Array{Float64}(undef,lines)
bcc_xa=Array{Float64}(undef,lines)
bcc_da=Array{Float64}(undef,lines)
bcc_dda=Array{Float64}(undef,lines)
bcc_db=Array{Float64}(undef,lines)
bcc_ddb=Array{Float64}(undef,lines)

file_id=open("bcc","r")
id=0
for i in eachline(file_id)
#	tempstring=readline(file_id)
	global id+=1
#	temp=split(i)[1]
	bcc_dens[id]=parse(Float64,split(i)[1])
#	println("$(temp) ||")
	bcc_xa[id]=parse(Float64,split(i)[2])
	bcc_da[id]=parse(Float64,split(i)[3])
	bcc_dda[id]=parse(Float64,split(i)[4])
	bcc_db[id]=parse(Float64,split(i)[5])
	bcc_ddb[id]=parse(Float64,split(i)[6])
end
close(file_id)

####################################333
file_id=open("fcc","r")
lines=0
for i in eachline(file_id)
	global lines+=1
end
close(file_id)
print("$(lines)")

fcc_dens=Array{Float64}(undef,lines)
fcc_xa=Array{Float64}(undef,lines)
fcc_da=Array{Float64}(undef,lines)
fcc_dda=Array{Float64}(undef,lines)
fcc_db=Array{Float64}(undef,lines)
fcc_ddb=Array{Float64}(undef,lines)

file_id=open("fcc","r")
id=0
for i in eachline(file_id)
#	tempstring=readline(file_id)
	global id+=1
#	temp=split(i)[1]
	fcc_dens[id]=parse(Float64,split(i)[1])
#	println("$(temp) ||")
	fcc_xa[id]=parse(Float64,split(i)[2])
	fcc_da[id]=parse(Float64,split(i)[3])
	fcc_dda[id]=parse(Float64,split(i)[4])
	fcc_db[id]=parse(Float64,split(i)[5])
	fcc_ddb[id]=parse(Float64,split(i)[6])
end
close(file_id)

#----get pure
i_a=0
i_b=0
for i=1:size(sc_dens,1)
	if(sc_xa[i]==1.0)
		global i_a+=1
	end
	if(sc_xa[i]==0.0)
		global i_b+=1
	end
end

sc_da=sc_da.*0.0109597281/0.3
sc_db=sc_db.*0.0109597281/0.3
bcc_da=bcc_da.*0.0109597281/0.3
bcc_db=bcc_db.*0.0109597281/0.3
fcc_da=fcc_da.*0.0109597281/0.3
fcc_db=fcc_db.*0.0109597281/0.3

sc_dda=sc_dda.*0.0109597281/0.3
sc_ddb=sc_dda.*0.0109597281/0.3
bcc_dda=bcc_dda.*0.0109597281/0.3
bcc_ddb=bcc_ddb.*0.0109597281/0.3
fcc_dda=fcc_dda.*0.0109597281/0.3
fcc_ddb=fcc_ddb.*0.0109597281/0.3

sc_pa_dens=Array{Float64}(undef,i_a)
sc_pa_d=Array{Float64}(undef,i_a)
sc_pa_dd=Array{Float64}(undef,i_a)

sc_pb_dens=Array{Float64}(undef,i_b)
sc_pb_d=Array{Float64}(undef,i_b)
sc_pb_dd=Array{Float64}(undef,i_b)

id_a=0
id_b=0
for i=1:size(sc_dens,1)
	if(sc_xa[i]==1.0)
		global id_a+=1
		sc_pa_dens[id_a]=sc_dens[i]
		sc_pa_d[id_a]=sc_da[i]
		sc_pa_dd[id_a]=sc_dda[i]
	end
	if(sc_xa[i]==0.0)
		global id_b+=1
		sc_pb_dens[id_b]=sc_dens[i]
		sc_pb_d[id_b]=sc_db[i]
		sc_pb_dd[id_b]=sc_ddb[i]
	end
end

file_id=open("sc_pa.our","w")
for i=1:i_a
	println(file_id,"$(sc_pa_dens[i])\t$(sc_pa_d[i])\t$(sc_pa_dd[i])")
end
close(file_id)

file_id=open("sc_pb.our","w")
for i=1:i_b
	println(file_id,"$(sc_pb_dens[i])\t$(sc_pb_d[i])\t$(sc_pb_dd[i])")
end
close(file_id)

#----get pure
i_a=0
i_b=0
for i=1:size(bcc_dens,1)
	if(bcc_xa[i]==1.0)
		global i_a+=1
	end
	if(bcc_xa[i]==0.0)
		global i_b+=1
	end
end
bcc_pa_dens=Array{Float64}(undef,i_a)
bcc_pa_d=Array{Float64}(undef,i_a)
bcc_pa_dd=Array{Float64}(undef,i_a)

bcc_pb_dens=Array{Float64}(undef,i_b)
bcc_pb_d=Array{Float64}(undef,i_b)
bcc_pb_dd=Array{Float64}(undef,i_b)

id_a=0
id_b=0
for i=1:size(bcc_dens,1)
	if(bcc_xa[i]==1.0)
		global id_a+=1
		bcc_pa_dens[id_a]=bcc_dens[i]
		bcc_pa_d[id_a]=bcc_da[i]
		bcc_pa_dd[id_a]=bcc_dda[i]
	end
	if(bcc_xa[i]==0.0)
		global id_b+=1
		bcc_pb_dens[id_b]=bcc_dens[i]
		bcc_pb_d[id_b]=bcc_db[i]
		bcc_pb_dd[id_b]=bcc_ddb[i]
	end
end

file_id=open("bcc_pa.our","w")
for i=1:i_a
	println(file_id,"$(bcc_pa_dens[i])\t$(bcc_pa_d[i])\t$(bcc_pa_dd[i])")
end
close(file_id)


file_id=open("bcc_pb.our","w")
for i=1:i_b
	println(file_id,"$(bcc_pb_dens[i])\t$(bcc_pb_d[i])\t$(bcc_pb_dd[i])")
end
close(file_id)

#----get pure
i_a=0
i_b=0
for i=1:size(fcc_dens,1)
	if(fcc_xa[i]==1.0)
		global i_a+=1
	end
	if(fcc_xa[i]==0.0)
		global i_b+=1
	end
end
fcc_pa_dens=Array{Float64}(undef,i_a)
fcc_pa_d=Array{Float64}(undef,i_a)
fcc_pa_dd=Array{Float64}(undef,i_a)

fcc_pb_dens=Array{Float64}(undef,i_b)
fcc_pb_d=Array{Float64}(undef,i_b)
fcc_pb_dd=Array{Float64}(undef,i_b)

id_a=0
id_b=0
for i=1:size(fcc_dens,1)
	if(fcc_xa[i]==1.0)
		global id_a+=1
		fcc_pa_dens[id_a]=fcc_dens[i]
		fcc_pa_d[id_a]=fcc_da[i]
		fcc_pa_dd[id_a]=fcc_dda[i]
	end
	if(fcc_xa[i]==0.0)
		global id_b+=1
		fcc_pb_dens[id_b]=fcc_dens[i]
		fcc_pb_d[id_b]=fcc_db[i]
		fcc_pb_dd[id_b]=fcc_ddb[i]
	end
end

file_id=open("fcc_pa.our","w")
for i=1:i_a
	println(file_id,"$(fcc_pa_dens[i])\t$(fcc_pa_d[i])\t$(fcc_pa_dd[i])")
end
close(file_id)


file_id=open("fcc_pb.our","w")
for i=1:i_b
	println(file_id,"$(fcc_pb_dens[i])\t$(fcc_pb_d[i])\t$(fcc_pb_dd[i])")
end
close(file_id)

#Read liquid data
liq_dens=Array{Float64}(undef,8)
liq_DA=Array{Float64}(undef,8)
liq_DDA=Array{Float64}(undef,8)
liq_DMA=Array{Float64}(undef,8)
liq_DDMA=Array{Float64}(undef,8)
liq_DB=Array{Float64}(undef,8)
liq_DDB=Array{Float64}(undef,8)
liq_DMB=Array{Float64}(undef,8)
liq_DDMB=Array{Float64}(undef,8)

file_id=open("liq.data","r")
id=0
for i=eachline(file_id)
	global id+=1
#	temp=split(i)[1]
	liq_dens[id]=parse(Float64,split(i)[1])
	liq_DA[id]=parse(Float64,split(i)[6])
	liq_DDA[id]=parse(Float64,split(i)[7])
	liq_DMA[id]=parse(Float64,split(i)[8])
	liq_DDMA[id]=parse(Float64,split(i)[9])

	liq_DB[id]=parse(Float64,split(i)[2])
	liq_DDB[id]=parse(Float64,split(i)[3])
	liq_DMB[id]=parse(Float64,split(i)[4])
	liq_DDMB[id]=parse(Float64,split(i)[5])
end
close(file_id)

liq_DA=liq_DA.*0.0109597281/0.3
liq_DB=liq_DB.*0.0109597281/0.3
liq_DDA=liq_DDA.*0.0109597281/0.3
liq_DDB=liq_DDB.*0.0109597281/0.3

liq_DMA=liq_DMA.*0.0109597281/0.3
liq_DMB=liq_DMB.*0.0109597281/0.3
liq_DDMA=liq_DDMA.*0.0109597281/0.3
liq_DDMB=liq_DDMB.*0.0109597281/0.3

file_id=open("comparison.out","w")
a1=bcc_pa_d./sc_pa_d
a2=fcc_pa_d./sc_pa_d
b1=bcc_pb_d./sc_pb_d
b2=fcc_pb_d./sc_pb_d

for i=1:i_b
	
	println(file_id,"$(sc_pa_dens[i])\t$(a1)\t$(a2)\t$(b1)\t$(b2)")
end
close(file_id)

#plot comparrisson
set_default_plot_size(16cm, 25cm) 

compbcc=layer(x=sc_pa_dens,y=a1, Geom.point, Theme(default_color="blue",key_position=:inside,point_size=1.5mm))
compfcc=layer(x=sc_pa_dens,y=a2, Geom.point, Theme(default_color="green",key_position=:inside,point_size=1.5mm))


c1=plot(compbcc,compfcc,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>/D<sup>*</sup><sub>SC</sub>"),Guide.title("A"),Guide.manual_color_key("", ["BCC", "FCC"], ["blue", "green"]))
#,Guide.manual_color_key("Тип решетки", ["BCC", "FCC"], ["blue", "green"])


compbcc=layer(x=sc_pa_dens,y=b1, Geom.point, Theme(default_color="blue", point_size=1.5mm))
compfcc=layer(x=sc_pa_dens,y=b2, Geom.point, Theme(default_color="green", point_size=1.5mm))

c2=plot(compbcc,compfcc,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>/D<sup>*</sup><sub>SC</sub>"),Guide.title("B"),Guide.manual_color_key("", ["BCC", "FCC"], ["blue", "green"]))

#,Guide.manual_color_key("Тип решетки", ["BCC", "FCC"], ["blue", "green"])

c3=vstack(c1,c2)
c3 |> PDF("comp.pdf")

#описание диффузии чистых веществ

@. pure(x, p) = p[1]*exp(x*p[2])+p[3]*exp(x*p[4])+p[5]

p0 = [0.5, 0.5, 0.5, 0.5, 0.5]

fit_pure_sc_a = curve_fit(pure, sc_pa_dens, sc_pa_d, p0)
println("sc a $(fit_pure_sc_a.param)")
fit_pure_sc_b = curve_fit(pure, sc_pb_dens, sc_pb_d, p0)
println("sc b $(fit_pure_sc_b.param)")

fit_pure_bcc_a = curve_fit(pure, bcc_pa_dens, bcc_pa_d, p0)
println("bcc a $(fit_pure_bcc_a.param)")
fit_pure_bcc_b = curve_fit(pure, bcc_pb_dens, bcc_pb_d, p0)
println("bcc b $(fit_pure_bcc_b.param)")

fit_pure_fcc_a = curve_fit(pure, fcc_pa_dens, fcc_pa_d, p0)
println("fcc a $(fit_pure_bcc_a.param)")
fit_pure_fcc_b = curve_fit(pure, fcc_pb_dens, fcc_pb_d, p0)
println("fcc b $(fit_pure_bcc_b.param)")

#plot results
#xdata = sort(iris[:SepalWidth])
#ydata = cumsum(xdata)
#line = layer(x=xdata, y=ydata, Geom.line, Theme(default_color="red"))
#bars = layer(iris, x=:SepalWidth, Geom.bar)
#plot(line, bars)

if(3==3)
xp=collect(0.05:0.01:0.8)
set_default_plot_size(16cm, 25cm)

pure_sc_a=layer(x=sc_pa_dens,y=sc_pa_d,ymin=(sc_pa_d.+sc_pa_dd), ymax=(sc_pa_d.-sc_pa_dd), Geom.point, Geom.errorbar, Theme(default_color="green"))
pure_sc_a_f=layer(x=xp,y=pure(xp,fit_pure_sc_a.param), Geom.line, Theme(default_color="green",line_style=[:solid],line_width=1.2mm))


pure_bcc_a=layer(x=bcc_pa_dens,y=bcc_pa_d,ymin=(bcc_pa_d.+bcc_pa_dd), ymax=(bcc_pa_d.-bcc_pa_dd), Geom.point, Geom.errorbar, Theme(default_color="red"))
pure_bcc_a_f=layer(x=xp,y=pure(xp,fit_pure_bcc_a.param), Geom.line, Theme(default_color="red",line_style=[:dash],line_width=1.2mm))

pure_fcc_a=layer(x=fcc_pa_dens,y=fcc_pa_d,ymin=(fcc_pa_d.+fcc_pa_dd), ymax=(fcc_pa_d.-fcc_pa_dd), Geom.point, Geom.errorbar, Theme(default_color="blue"))
pure_fcc_a_f=layer(x=xp,y=pure(xp,fit_pure_fcc_a.param), Geom.line, Theme(default_color="blue",line_style=[:dashdot],line_width=1.2mm))

liq_a=layer(x=liq_dens,y=liq_DA,ymin=(liq_DA.+liq_DDA), ymax=(liq_DA.-liq_DDA), Geom.point, Geom.errorbar, Theme(default_color="red"))

p1=plot(pure_sc_a,pure_sc_a_f,pure_bcc_a,pure_bcc_a_f,pure_fcc_a,pure_fcc_a_f,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>"),Guide.title("A"))

pp1=plot(pure_sc_a,pure_sc_a_f,pure_bcc_a,pure_bcc_a_f,pure_fcc_a,pure_fcc_a_f,liq_a,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>"),Guide.title("A"))

#b
pure_sc_b=layer(x=sc_pb_dens,y=sc_pb_d,ymin=(sc_pb_d.+sc_pb_dd), ymax=(sc_pb_d.-sc_pb_dd), Geom.point, Geom.errorbar, Theme(default_color="green"))
pure_sc_b_f=layer(x=xp,y=pure(xp,fit_pure_sc_b.param), Geom.line, Theme(default_color="green",line_style=[:solid],line_width=1.2mm))

pure_bcc_b=layer(x=bcc_pb_dens,y=bcc_pb_d,ymin=(bcc_pb_d.+bcc_pb_dd), ymax=(bcc_pb_d.-bcc_pb_dd), Geom.point, Geom.errorbar, Theme(default_color="red"))
pure_bcc_b_f=layer(x=xp,y=pure(xp,fit_pure_bcc_b.param), Geom.line, Theme(default_color="red",line_style=[:dash],line_width=1.2mm))

pure_fcc_b=layer(x=fcc_pb_dens,y=fcc_pb_d,ymin=(fcc_pb_d.+fcc_pb_dd), ymax=(fcc_pb_d.-fcc_pb_dd), Geom.point, Geom.errorbar, Theme(default_color="blue"))
pure_fcc_b_f=layer(x=xp,y=pure(xp,fit_pure_fcc_b.param), Geom.line, Theme(default_color="blue",line_style=[:dashdot],line_width=1.2mm))

liq_b=layer(x=liq_dens,y=liq_DB,ymin=(liq_DB.+liq_DDB), ymax=(liq_DB.-liq_DDB), Geom.point, Geom.errorbar, Theme(default_color="red"))

#p1=plot(pure_sc_a,pure_sc_a_f,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>"))
p2=plot(pure_sc_b,pure_sc_b_f,pure_bcc_b,pure_bcc_b_f,pure_fcc_b,pure_fcc_b_f,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>"),Guide.title("B"))

pp2=plot(pure_sc_b,pure_sc_b_f,pure_bcc_b,pure_bcc_b_f,pure_fcc_b,pure_fcc_b_f, liq_b ,Guide.xlabel("n<sup>*</sup>"), Guide.ylabel("D<sup>*</sup>"),Guide.title("B"))

p3=vstack(p1,p2)
pp3=vstack(pp1,pp2)
p3 |> PDF("PureAB.pdf")
pp3 |> PDF("PureABv2.pdf")

end


#scatter(sc_pa_dens,sc_pa_d,label="simple cubic")
#scatter!(bcc_pa_dens,bcc_pa_d,label="bcc")
#scatter!(fcc_pa_dens,fcc_pa_d,label="fcc")
#xlabel!("$n^*$")
#ylabel!("$D^*$")
#savefig("pure.pdf")

#####################------------Mixture

#@. MixDscA(x,p)=pure(x[1],fit_pure_sc_a.param)+p*x[1]*(1-x[2])

function MixDscA(x,p)
	#println("   x31  $(x[3,1])   x31 $(floor(Int,x[3,1]))")
	
	out=Array{Float64}(undef,floor(Int,x[3,1]))
	for i=1:floor(Int,x[3,1])
		out[i]=abs(x[5,i] -pure(x[1,i],x[4,:])+p[1]*x[1,i]*(1.0-x[2,i]))/x[5,i]
	end
	return vec(out)
end

function MixD(x,p)
	#println("   x31  $(x[3,1])   x31 $(floor(Int,x[3,1]))")
	
	out=Array{Float64}(undef,floor(Int,x[3,1]))
	for i=1:floor(Int,x[3,1])
		out[i]=pure(x[1,i],x[4,:])+p[1]*x[1,i]*(1.0-x[2,i])
	end
	return vec(out)
end

##CUBIC A
id=0 
for i=1:size(sc_da,1)
	if sc_da[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if sc_da[i]>0
		global id+=1
		xin[1,id]=sc_dens[i]
		xin[2,id]=sc_xa[i]
		xin[5,id]=sc_da[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_sc_a.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_sc_a = curve_fit(MixDscA, xin, yin, p0)
println("sc a $(fit_m_sc_a.param)")



test=MixD(xin,fit_m_sc_a.param)
test=MixD(xin,[0.00656768])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("sc_a.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

##CUBIC B
id=0 
for i=1:size(sc_db,1)
	if sc_db[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if sc_db[i]>0
		global id+=1
		xin[1,id]=sc_dens[i]
		xin[2,id]=1.0-sc_xa[i]
		xin[5,id]=sc_db[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_sc_b.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_sc_b = curve_fit(MixDscA, xin, yin, p0)
println("sc b $(fit_m_sc_b.param)")

test=MixD(xin,fit_m_sc_b.param)
test=MixD(xin,[-0.00656768])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("sc_b.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

pa=(fit_m_sc_a.param.-fit_m_sc_b.param)./4
pb=-(fit_m_sc_a.param.-fit_m_sc_b.param)./4
println("pa  $(pa) pb $(pb)")

##BCC A
id=0 
for i=1:size(bcc_da,1)
	if bcc_da[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if bcc_da[i]>0
		global id+=1
		xin[1,id]=bcc_dens[i]
		xin[2,id]=bcc_xa[i]
		xin[5,id]=bcc_da[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_bcc_a.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_bcc_a = curve_fit(MixDscA, xin, yin, p0)
println("bcc a $(fit_m_bcc_a.param)")



test=MixD(xin,fit_m_bcc_a.param)
test=MixD(xin,[0.000407613])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("bcc_a.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

##CUBIC B
id=0 
for i=1:size(bcc_db,1)
	if bcc_db[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if bcc_db[i]>0
		global id+=1
		xin[1,id]=bcc_dens[i]
		xin[2,id]=1.0-bcc_xa[i]
		xin[5,id]=bcc_db[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_bcc_b.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_bcc_b = curve_fit(MixDscA, xin, yin, p0)
println("bcc b $(fit_m_bcc_b.param)")

test=MixD(xin,fit_m_bcc_b.param)
test=MixD(xin,[-0.000407613])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("bcc_b.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

pa=(fit_m_bcc_a.param.-fit_m_bcc_b.param)./2
pb=-(fit_m_bcc_a.param.-fit_m_bcc_b.param)./2
println("pa  $(pa) pb $(pb)")

##FCC A
id=0 
for i=1:size(fcc_da,1)
	if fcc_da[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if fcc_da[i]>0
		global id+=1
		xin[1,id]=fcc_dens[i]
		xin[2,id]=fcc_xa[i]
		xin[5,id]=fcc_da[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_fcc_a.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_fcc_a = curve_fit(MixDscA, xin, yin, p0)
println("fcc a $(fit_m_fcc_a.param)")



test=MixD(xin,fit_m_fcc_a.param)
test=MixD(xin,[0.00180914])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("fcc_a.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

##CUBIC B
id=0 
for i=1:size(fcc_db,1)
	if fcc_db[i]>0
		global id+=1;
	end
end
dim=id
xin=Array{Float64}(undef,5,dim)
yin=Array{Float64}(undef,dim)
id=0
for i=1:lines
	if fcc_db[i]>0
		global id+=1
		xin[1,id]=fcc_dens[i]
		xin[2,id]=1.0-fcc_xa[i]
		xin[5,id]=fcc_db[i]
		yin[id]=0	#sc_da[i]
	end
end
xin[3,1]=dim
for i=1:5
	xin[4,i]=fit_pure_fcc_b.param[i]
end

#println(" $(xin)")

p0=[0.02]
test=MixDscA(xin,p0)
#println("   test $(test)")
fit_m_fcc_b = curve_fit(MixDscA, xin, yin, p0)
println("fcc b $(fit_m_fcc_b.param)")

test=MixD(xin,fit_m_fcc_b.param)
test=MixD(xin,[-0.00180914])
err1=test.-xin[5,:]
err2=(test.-xin[5,:])./xin[5,:].*100


file_id=open("fcc_b.test","w")
for i=1:dim
	println(file_id,"$(xin[1,i])\t$(xin[2,i])\t$(xin[5,i])\t$(test[i])\t$(err2[i])")
end
close(file_id)

pa=(fit_m_fcc_a.param.-fit_m_fcc_b.param)./2
pb=-(fit_m_fcc_a.param.-fit_m_fcc_b.param)./2
println("pa  $(pa) pb $(pb)")

