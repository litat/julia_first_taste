using DataFrames, GLM, Gadfly

cd("/Users/lita/.jl/first_taste/")
glossy = readtable("./glossy_all.csv")

pvalue(lm1) = [ccdf(FDist(1,df_residual(lm1.model)),abs2(fval)) for fval in coef(lm1)./stderr(lm1)][2]

for treat = unique(glossy[:treatments])
	glossytreat = glossy[glossy[:treatments] .== treat, :]
	lm1 = lm(glossy~days, glossytreat)
	if pvalue(lm1) < 0.05
		println("\nGlossy differ in days by $treat")
		println(lm1)
	end
end

blank = glossy[(glossy[:treatments] .== "Blank"), :]
btq = glossy[(glossy[:treatments] .== "Blank") | (glossy[:treatments] .== "T1130") | (glossy[:treatments] .== "Quercetin"), :]
ncg = glossy[(glossy[:treatments] .== "NC") | (glossy[:treatments] .== "NCT") | (glossy[:treatments] .== "NC+T1130") | (glossy[:treatments] .== "NC+Quercetin"), :]
pug = glossy[(glossy[:treatments] .== "PU") | (glossy[:treatments] .== "PUT") | (glossy[:treatments] .== "PU+T1130") | (glossy[:treatments] .== "PU+Quercetin"), :]

for group = [[btq] [ncg] [pug]], day = unique(glossy[:days])
	d = group[group[:days] .== day, :]
	pool!(d)
	lm1d = lm(glossy~treatments, d)
	println("\nGlossy differ in treatments by $day")
	println(lm1d)
end

ncpu = glossy[(glossy[:treatments] .== "NC") | (glossy[:treatments] .== "PU"), :]
nctput = glossy[(glossy[:treatments] .== "NCT") | (glossy[:treatments] .== "PUT"), :]
ncatpuat = glossy[(glossy[:treatments] .== "NC+T1130") | (glossy[:treatments] .== "PU+T1130"), :]
ncqpuq = glossy[(glossy[:treatments] .== "NC+Quercetin") | (glossy[:treatments] .== "PU+Quercetin"), :]

for group = [[ncpu] [nctput] [ncatpuat] [ncqpuq]], day = unique(glossy[:days])
	d = group[group[:days] .== day, :]
	pool!(d)
	lm1dt = lm(glossy~treatments, d)
	println("\nGlossy differ in treatments by $day")
	println(lm1dt)
end

set_default_plot_size(30cm, 18cm)

plot_glossy = plot(glossy, xgroup="treatments", x="days", y="glossy", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_glossy.svg", 30cm, 18cm), plot_glossy)

plot_glossy_blank = plot(blank, xgroup="orientation", x="glossy", color="days", Geom.subplot_grid(Geom.density))
draw(SVG("plot_glossy_blank.svg", 30cm, 18cm), plot_glossy_blank)

set_default_plot_size(10cm, 18cm)
plot_glossy_btq = plot(btq, xgroup="treatments", x="days", y="glossy", Scale.x_discrete, Geom.subplot_grid(Geom.boxplot))
draw(SVG("plot_glossy_btq.svg", 10cm, 18cm), plot_glossy_btq)


set_default_plot_size(12.5cm, 18cm)

plot_glossy_ncg = plot(ncg, xgroup="days", x="treatments", y="glossy", Scale.x_discrete, Geom.subplot_grid(Geom.boxplot))
draw(SVG("plot_glossy_ncg.svg", 10cm, 18cm), plot_glossy_ncg)

plot_glossy_pug = plot(pug, xgroup="days", x="treatments", y="glossy", Scale.x_discrete, Geom.subplot_grid(Geom.boxplot))
draw(SVG("plot_glossy_pug.svg", 10cm, 18cm), plot_glossy_pug)

