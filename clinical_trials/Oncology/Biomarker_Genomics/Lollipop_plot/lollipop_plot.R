library(ggplot2); library(ggrepel); library(dplyr)

set.seed(1)
n <- 38
tp53 <- data.frame(
  USUBJID     = paste0("PT",sprintf("%03d",1:n)),
  PROTEIN_POS = sample(1:1100, n),
  VAF         = round(runif(n, 0.08, 0.78), 2),
  VARIANT_CLASS = sample(c("Missense_Mutation","Nonsense_Mutation",
                            "Frame_Shift_Del","Splice_Site"), n, TRUE,
                          prob=c(.55,.2,.15,.1)),
  CLONAL      = sample(c("Y","N"), n, TRUE, prob=c(.65,.35)),
  TUMORTYPE   = sample(c("NSCLC","BRCA","HCC","CRC","PDAC"), n, TRUE),
  HGVSP       = paste0("p.",sample(c("R175H","G245S","R248W","R273H",
                                      "P309S","V143A"), n, TRUE))
)

VC_COLORS <- c(Missense_Mutation="#2c5f8a", Nonsense_Mutation="#c0392b",
               Frame_Shift_Del="#e67e22", Splice_Site="#8e44ad")
TT_COLORS <- c(NSCLC="#2c5f8a", BRCA="#c0392b", HCC="#27ae60",
               CRC="#e67e22",   PDAC="#8e44ad")

lbl <- tp53 %>% filter(VAF >= 0.55) %>%
  mutate(label = sub("p\\.","",HGVSP))

p <- ggplot(tp53, aes(x=PROTEIN_POS, y=VAF)) +
  annotate("rect", xmin=100, xmax=1100, ymin=-0.025, ymax=0, fill="#2c5f8a") +
  annotate("text", x=600, y=-0.013, label="DNA-Binding Domain",
           color="white", size=3, fontface="bold") +
  geom_segment(aes(xend=PROTEIN_POS, yend=0, color=VARIANT_CLASS),
               linewidth=0.8, alpha=0.7) +
  geom_point(data=tp53 %>% filter(CLONAL=="Y"),
             aes(color=VARIANT_CLASS), size=3.5) +
  geom_point(data=tp53 %>% filter(CLONAL=="N"),
             aes(color=VARIANT_CLASS), size=3.5, shape=1, stroke=1.2) +
  geom_text_repel(data=lbl, aes(label=label), size=2.8,
                  box.padding=0.4, segment.color="#bbbbbb") +
  scale_color_manual(values=c(VC_COLORS, TT_COLORS), na.value="#888") +
  scale_x_continuous(name="Amino Acid Position", breaks=seq(0,1100,200)) +
  scale_y_continuous(name="VAF", limits=c(-0.05, 0.85)) +
  theme_classic(base_size=10) +
  labs(title="Lollipop Plot — TP53 Somatic Mutations",
       subtitle="n = 38 mutations · filled = clonal, open = subclonal",
       color="Variant Class")

ggsave("plots/01_lollipop.png", p, width=16, height=7.5, dpi=180, bg="white")
cat("01 done\n")
