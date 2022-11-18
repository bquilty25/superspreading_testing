a <- ggplot(data=data.frame(x=0:20))+
  stat_function(aes(x=x),fun=dnbinom,args=list(mu=2.5,size=100),geom="bar",n=21,fill=bi_col_pal[1])+
  labs(x="Secondary infections (R)",y="Density",title="R: 2.5, k: 100")+
  plotting_theme+
  ylim(0,0.8)


b <- ggplot(data=data.frame(x=0:20))+
  stat_function(aes(x=x),fun=dnbinom,args=list(mu=2.5,size=0.1),geom="bar",n=21,fill=bi_col_pal[2])+
  labs(x="Secondary infections (R)",y="Density",title="R: 2.5, k: 0.1")+
  plotting_theme+
  ylim(0,0.8)

a
b

ggsave(a,filename="results/high_k.png",width=100,height=100,units="mm",dpi=600)
ggsave(b,filename="results/low_k.png",width=100,height=100,units="mm",dpi=600)

