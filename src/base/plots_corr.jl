using TidierFiles
using TidierPlots

ads = read_csv("https://raw.githubusercontent.com/ivaquero/blog-bio/main/analysis/data/advertising.csv")

p1 = ggplot(ads, @aes(x = TV, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5);

p2 = ggplot(ads, @aes(x = radio, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5);

p3 = ggplot(ads, @aes(x = newspaper, y = sales)) +
     geom_point(color=:red) +
     geom_smooth(method=:lm, linewidth=3, alpha=0.5);

p1 + p2 + p3 + plot_layout(ncol=3, nrow=1, widths=[2, 2, 2], heights=[1])
