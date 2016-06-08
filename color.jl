using DataFrames, GLM, Gadfly

cd("/Users/lita/.jl/first_taste/")
color = readtable("color_all.csv")

L(Y) = 116*cbrt(Y/100)-16
a(X, Y) = 500(cbrt(X/94.81)-cbrt(Y/100))
b(Y, Z) = 200(cbrt(Y/100)-cbrt(Z/107.34))
YI(X, Y, Z) = 100(1.28X-1.06Z)./Y
C(a, b) = sqrt(a.^2+b.^2)
H(a, b) = atan(b./a)

color[:L] = L(color[:Y0])
color[:a] = a(color[:X0], color[:Y0])
color[:b] = b(color[:Y0], color[:Z0])
color[:YI] = YI(color[:X0], color[:Y0], color[:Z0])
color[:C] = C(color[:a], color[:b])
color[:H] = H(color[:a], color[:b])

pvalue(lm1) = [ccdf(FDist(1,df_residual(lm1.model)),abs2(fval)) for fval in coef(lm1)./stderr(lm1)][2]

# YI
for treat = unique(color[:treatments])
	colortreat = color[color[:treatments] .== treat, :]
	lm1 = lm(YI~days, colortreat)
	# if pvalue(lm1) < 0.05
		println("\nYI differ in days by $treat")
		println(lm1)
	# end
end

# C
for treat = unique(color[:treatments])
	colortreat = color[color[:treatments] .== treat, :]
	lm1 = lm(C~days, colortreat)
	# if pvalue(lm1) < 0.05
		println("\nC differ in days by $treat")
		println(lm1)
	# end
end

# H
for treat = unique(color[:treatments])
	colortreat = color[color[:treatments] .== treat, :]
	lm1 = lm(H~days, colortreat)
	# if pvalue(lm1) < 0.05
		println("\nH differ in days by $treat")
		println(lm1)
	# end
end



blank = color[(color[:treatments] .== "Blank"), :]
btq = color[(color[:treatments] .== "Blank") | (color[:treatments] .== "T1130") | (color[:treatments] .== "Quercetin"), :]
ncg = color[(color[:treatments] .== "NC") | (color[:treatments] .== "NCT") | (color[:treatments] .== "NC+T1130") | (color[:treatments] .== "NC+Quercetin"), :]
pug = color[(color[:treatments] .== "PU") | (color[:treatments] .== "PUT") | (color[:treatments] .== "PU+T1130") | (color[:treatments] .== "PU+Quercetin"), :]

for group = [[btq] [ncg] [pug]], day = unique(color[:days])
	d = group[group[:days] .== day, :]
	pool!(d)
	lm1d = lm(YI~treatments, d)
	println("\nYI differ in treatments by $day")
	println(lm1d)
end

for group = [[btq] [ncg] [pug]], day = unique(color[:days])
	d = group[group[:days] .== day, :]
	pool!(d)
	lm1d = lm(C~treatments, d)
	println("\nC differ in treatments by $day")
	println(lm1d)
end

for group = [[btq] [ncg] [pug]], day = unique(color[:days])
	d = group[group[:days] .== day, :]
	pool!(d)
	lm1d = lm(H~treatments, d)
	println("\nH differ in treatments by $day")
	println(lm1d)
end



set_default_plot_size(30cm, 18cm)
plot_color_L = plot(color, xgroup="treatments", x="days", y="L", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_L.svg", 30cm, 18cm), plot_color_L)

plot_color_a = plot(color, xgroup="treatments", x="days", y="a", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_a.svg", 30cm, 18cm), plot_color_a)

plot_color_b = plot(color, xgroup="treatments", x="days", y="b", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_b.svg", 30cm, 18cm), plot_color_b)

plot_color_YI = plot(color, xgroup="treatments", x="days", y="YI", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_YI.svg", 30cm, 18cm), plot_color_YI)

plot_color_C = plot(color, xgroup="treatments", x="days", y="C", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_C.svg", 30cm, 18cm), plot_color_C)

plot_color_H = plot(color, xgroup="treatments", x="days", y="H", Geom.subplot_grid(Geom.boxplot), Scale.x_discrete)
draw(SVG("plot_color_H.svg", 30cm, 18cm), plot_color_H)

set_default_plot_size(30cm, 18cm)

plot_color_days = plot(color, xgroup="days", x="treatments", y="YI", Scale.x_discrete, Geom.subplot_grid(Geom.boxplot))
draw(SVG("plot_color_days.svg", 30cm, 18cm), plot_color_days)

set_default_plot_size(30cm, 100cm)
plot_color_ab = plot(color, xgroup="days", ygroup="treatments", x="a", y="b", color="days", Geom.subplot_grid(Geom.point))
draw(SVG("plot_color_ab.svg", 35cm, 100cm), plot_color_ab)

draw(SVG("plot_color_ba.svg", 90cm, 36cm), plot(color, xgroup="treatments", ygroup="days", x="a", y="b", color="days", Geom.subplot_grid(Geom.point)))
draw(SVG("plot_color_treatments.svg", 15cm, 140cm), plot(color, ygroup="treatments", x="a", y="b", color="days", Geom.subplot_grid(Geom.point)))
