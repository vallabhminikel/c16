# STARTUP #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

# DEPENDENCIES ####

tell_user('Loading required packages...')

options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/c16')
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(survival)
  library(janitor)
  library(magick)
  library(smoother)
  select = dplyr::select
  summarize = dplyr::summarize
})


# OUTPUT STREAMS #### 

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/data-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}



# FUNCTIONS ####

tell_user('done.\nDefining functions...')


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

rank_uniq = function(x) match(x, sort(unique(x)))

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}

lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

repeach = function(vec, each) { # because R's default rep() only takes 1 value for "each"
  repped = c()
  for (i in 1:length(vec)) {
    repped = c(repped, rep(vec[i],each=each[i]))
  }
  repped
}

rbind_files = function(path, grepstring, keep_cols=NULL) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    if (!is.null(keep_cols)) {
      this_tbl = this_tbl[,keep_cols]
    }
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}


## Spearman-Karber functions from https://www.cureffi.org/2015/09/20/the-math-behind-spearman-karber-analysis/ ####

smooth_p = function(p) {
  if (length(p)<2) {
    return (p)
  }
  while (any(cummin(p) != p)) { # while p is non-monotonic
    for (i in 2:length(p)) { # starting from the second dilution and going up,
      if (p[i] > p[i-1]) { # if the current one has higher p than the last one...
        # select all dilutions below the current one that have lower p...
        indices_to_average = c(which(1:length(p) < i & p < p[i]), i)
        # and average them all together with the current value.
        p[indices_to_average] = mean(p[indices_to_average])
      }
    }
  }
  return (p)
}


spearman_karber = function(x, n, d) {
  if (all(x==0)) { # samples negative at a single tested dilution = LLQ
    estimate = -d[1]
    se = 0
  } else {
    if (length(x)==1) {
      dd = 1
    } else {
      dd = min(abs(diff(d))) # find dd, the difference between dilutions
      if (any(abs(diff(d)) - dd > .001)) { # .001 provides some numeric tolerance
        stop("x are not uniformly distributed\n")
      }
    }
    # if a non-zero positive is observed at a dilution higher than where zeroes are observed, call this noise and remove it
    for (j in 2:length(x)) {
      if (x[j-1]==0 & x[j] > 0) {
        x[j] = 0
      }
    }
    nd = length(x) # number of dilutions
    p = smooth_p(x/n) # smooth the proportion positive to be non-increasing
    if (all(p < 1)) { # if all positive responses never observed
      i = 1
      # first, imagine that it would have been all positive at one higher dilution
      sd50 = (d[i]+dd) - dd*(1/2 + sum(p[i:nd])) # calculate log10(sd50)
      # but never report a value lower than LLQ
      sd50 = pmin(sd50, d[i])
      estimate = -sd50 # estimate of seeds in original sample is opposite of the sd50
      se = sqrt(dd^2 * sum(p*(1-p)/(n-1))) # calculate standard error
    } else {
      i = max(which(p==1)) # find the weakest dilution with all positive responses
      sd50 = d[i] - dd*(1/2 + sum(p[(i+1):nd])) # calculate log10(sd50)
      estimate = -sd50 # estimate of seeds in original sample is opposite of the sd50
      se = sqrt(dd^2 * sum(p*(1-p)/(n-1))) # calculate standard error
    }
  }
  return_value = cbind(estimate=estimate, se=se)
  return (return_value)
}

## survival wraparounds ####

make_sf = function(survdata) {
  sf = survfit(Surv(dpi, acm) ~ cohort, data=survdata)
  sf$cohort = gsub('cohort=','',names(sf$strata))
  sf$color = survdata$color[match(sf$cohort, survdata$cohort)]
  sf$disp = survdata$disp[match(sf$cohort, survdata$cohort)]
  sf$lty = survdata$lty[match(sf$cohort, survdata$cohort)]
  return(sf)
}

make_sv = function(survdata) {
  sv = survdiff(Surv(dpi, acm) ~ cohort, data=survdata)
  sv$cohort = gsub('cohort=','',names(sv$strata))
  sv$color = survdata$color[match(sv$cohort, survdata$cohort)]
  sv$disp = survdata$disp[match(sv$cohort, survdata$cohort)]
  sv$lty = survdata$lty[match(sv$cohort, survdata$cohort)]
  return(sv)
}

sv_pval = function(sv) {
  return (pchisq(q=sv$chisq, df=length(sv$n)-1, lower.tail=F))
}

sf_table = function(sf) {
  # per-cohort N, events, and median survival table for the supplement, from a survfit object
  tab = summary(sf)$table
  if (is.null(dim(tab))) tab = matrix(tab, nrow=1, dimnames=list(names(sf$strata), names(tab)))
  as_tibble(tab, rownames='stratum') %>%
    clean_names() %>%
    mutate(cohort = gsub('cohort=','',stratum)) %>%
    select(-stratum) %>%
    left_join(tibble(cohort=sf$cohort, disp=sf$disp) %>% distinct(), by='cohort') %>%
    select(cohort, disp, records, n_start, events, median, x0_95lcl, x0_95ucl)
}


plot_sf = function(sf, 
                   anno_dpi=NULL, 
                   xlims=c(0,500), 
                   xats = seq(0,500,50),
                   xbigs = seq(0,500,100),
                   ylims = c(0,1.05),
                   yats = seq(0,1,.25),
                   show_legend=T) {
  
  if(!('lty' %in% names(sf))) sf$lty=1
  
  plot(sf, col=sf$color, lty=sf$lty, lwd=2, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  censored = sf$n.censor > 0 & sf$time < 500
  points(sf$time[censored], sf$surv[censored], pch=1, lwd=1, col=repeach(sf$color,each=sf$strata)[censored])
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.6, text='days post-inoculation', cex=0.8)
  axis(side=2, at=yats, labels=NA, tck=-0.05)
  axis(side=2, at=yats, labels=percent(yats), lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='survival', cex=0.8)
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    points(x=c(anno_dpi), y=rep(max(ylims),length(anno_dpi)), pch=25, bg='#000000')
  }
  if(show_legend) {
    leg = tibble(disp=sf$disp, color=sf$color, lty=sf$lty) %>% 
      distinct(disp, color, lty)
    legend(x=max(xlims), y=1, bty='n', cex=.8, lwd=2, legend=leg$disp,
           col=leg$color, lty=leg$lty,
           text.col=leg$color)
  }
  par(xpd=F)
}

norm_wts = function(weights_raw, baseline_dpi, main_tbl, meta_tbl) {
  weights_raw %>%
    filter(wtdpi %in% baseline_dpi) %>%
    select(animal, baseline=weight) -> baselines
  weights_raw %>%
    inner_join(baselines, by='animal') %>%
    mutate(rel=weight/baseline-1) %>%
    inner_join(select(main_tbl, animal, cohort), by='animal') %>%
    inner_join(meta_tbl, by=c('cohort'))  -> rel_weights
  rel_weights %>%
    group_by(cohort, disp, color, lty, wtdpi) %>%
    summarize(.groups='keep',
              mean = mean(rel),
              sd = sd(rel),
              n = n(),
              l95 = mean(rel) - 1.96*sd(rel)/sqrt(n()),
              u95 = mean(rel) + 1.96*sd(rel)/sqrt(n())) %>%
    filter(n > 1) -> rel_weights_smry
  return (list(rel_weights, rel_weights_smry))
}



plot_wts = function(normed_weights_list, 
                    anno_dpi=NULL, 
                    xlims=c(0,500), 
                    xats = seq(0,500,50),
                    xbigs = seq(0,500,100),
                    ylims = c(-0.25,0.4),
                    yats = seq(-0.3,0.5,0.1),
                    show_legend=F) {
  
  rel_weights_smry = normed_weights_list[[2]]
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.6, text='days post-inoculation', cex=0.8)
  axis(side=2, at=yats, labels=NA, tck=-0.05)
  axis(side=2, at=yats, labels=percent(yats), lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='weight (% change)', cex=0.8)
  abline(h=0)
  abline(h=-.2, col='red', lty=3)
  mtext(side=4, at=-.2, las=2, col='red', text='endpoint', cex=.8)
  for (this_cohort in unique(rel_weights_smry$cohort)) {
    subs = rel_weights_smry %>% filter(cohort==this_cohort) 
    points(subs$wtdpi, subs$mean, type='l', lwd=2, col=subs$color)
    polygon(x=c(subs$wtdpi,rev(subs$wtdpi)),y=c(subs$l95,rev(subs$u95)),col=alpha(subs$color,ci_alpha),border=NA)
  }
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    points(x=c(anno_dpi), y=rep(max(ylims),length(anno_dpi)), pch=25, bg='#000000')
  }
  if(show_legend) {
    leg = tibble(disp=rel_weights_smry$disp, color=rel_weights_smry$color, lty=rel_weights_smry$lty) %>% distinct(disp, color, lty)
    legend(x=max(xlims), y=1, bty='n', cex=.8, lwd=2, legend=leg$disp,
           col=leg$color, lty=leg$lty,
           text.col=leg$color)
  }
  par(xpd=F)
}

norm_nests = function(nests_raw, main_tbl) {
  main_tbl %>%
    select(cage, cohort, color, disp, lty) %>%
    group_by(cage) %>%
    slice(1) %>%
    ungroup() -> cage_corresp

  cage_corresp %>%
    inner_join(nests_raw, by='cage') %>%
    rename(mean=comb) -> nests_anno
  
  return (nests_anno)
}

plot_nests = function(nests_anno, 
                      anno_dpi=NULL, 
                      xlims=c(0,500), 
                      xats = seq(0,500,50),
                      xbigs = seq(0,500,100),
                      ylims = c(0,2.05),
                      yats = seq(0,2,0.5),
                      show_legend=F)  {
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.6, text='days post-inoculation', cex=0.8)
  axis(side=2, at=yats, labels=NA, tck=-0.05)
  axis(side=2, at=yats, labels=yats, lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='nest score', cex=0.8)
  
  # thin line for each cage
  for (current_cage in unique(nests_anno$cage)) {
    subs = nests_anno %>% filter(cage==current_cage) 
    points(subs$nestdpi, subs$mean, type='l', col=subs$color[1],lwd=0.25)
  }
  
  # trend line for each group
  for (this_cohort in unique(nests_anno$cohort)) {
    subs = nests_anno %>% 
      filter(cohort==this_cohort) %>% 
      group_by(cohort,nestdpi,color) %>% 
      summarise(.groups='keep',
                n=n(),
                mean=mean(mean)) %>%
      ungroup()
    gaussian_smoothed = smth(subs$mean, method = "gaussian", window = 10, tails = TRUE)
    points(subs$nestdpi, gaussian_smoothed, type = 'l', col = subs$color[1], lwd = 2)
    
  }
  
  nests_anno %>% distinct(disp, color) -> leg
  
  par(xpd=T)
  if (!is.null(anno_dpi)) {
    points(x=c(anno_dpi), y=rep(max(ylims),length(anno_dpi)), pch=25, bg='#000000')
  }
  if(show_legend) {
    legend(x=max(xlims), y=1, bty='n', cex=.8, lwd=2, legend=leg$disp,
           col=leg$color, lty=leg$lty,
           text.col=leg$color)
  }
  par(xpd=F)
  
  
}


# DATA ####

tell_user('done.\nReading in data...')

## Analytic datasets ####

activity = read_tsv('data/activity.tsv',col_types=cols())
alny_data = read_tsv('data/alny_data.tsv',col_types=cols())
blind = read_tsv('data/blind.tsv',col_types=cols())
elisa = read_tsv('data/elisa_6d11.tsv',col_types=cols())
main = read_tsv('data/main.tsv',col_types=cols())
meta = read_tsv('data/meta.tsv',col_types=cols())
nests = read_tsv('data/nests.tsv',col_types=cols())
nfl = read_tsv('data/nfl.tsv',col_types=cols())
qpcr = read_tsv('data/qpcr.tsv', col_types=cols())
rtquics = read_tsv('data/rtquics.tsv',col_types=cols())
weights = read_tsv('data/weights.tsv',col_types=cols())

opt_6d11 = read_tsv('supplemental/6D11_ELISA_optimization.tsv', col_types=cols())

## Process RT-QuIC titers ####
bh_volume_ml = 0.002 # 2 uL seed. adjust for this in titer calculation to get units of SD50 per gram at end
rtquics %>%
  group_by(sample) %>%
  arrange(sample, desc(dilution)) %>%
  summarize(.groups='keep',
            sk = spearman_karber(x=positive_replicates, n=total_replicates, d=log10(dilution))) %>% 
  ungroup() %>%
  mutate(titer = sk[,1] + -log10(bh_volume_ml), se=sk[,2]) %>%
  select(-sk) -> rtquic_titers

## Merge survival data ####
survdata = main %>% 
  select(-cage) %>%
  inner_join(blind, by=c('animal','expt')) %>%
  inner_join(meta, by='cohort') %>%
  filter(!exclude) %>%
  select(-exclude)

# DISPLAY ITEMS ####

## Supplemental data ####

### age distribution ####
main %>%
  mutate(age_at_inoc = as.integer(inoculation_date - dob),
         dpi_at_baseline = 112- as.integer(inoculation_date - dob)) %>%
  summarize(.groups='keep',
            min_weeks_at_inoc=min(age_at_inoc, na.rm=T)/7,
            max_weeks_at_inoc=max(age_at_inoc, na.rm=T)/7,
            min_dpi_at_baseline=min(dpi_at_baseline, na.rm=T),
            max_dpi_at_baseline=max(dpi_at_baseline, na.rm=T)) %>%
  ungroup() -> age_range 
write_supp_table(age_range, 'Range of ages at inoculation and dpi at baseline for Methods section.')

### sex distribution ####
main %>%
  group_by(expt) %>%
  summarize(.groups='keep',
            n_m = sum(sex=='M'), 
            n_f = sum(sex=='F')
            ) %>%
  ungroup() -> sex_distribution 
write_supp_table(sex_distribution, 'Number of animals by sex in each experiment.')



## Figure S1.  #### 
tell_user('done.\nCreating Figure S1...')
if ('figure-s1'=='figure-s1') {
  
  resx=300
  png('display_items/figure-s1.png',width=6.5*resx,height=4.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,2), nrow=1, byrow=T)
  layout(layout_matrix)
  
  panel = 1

  ### A. Endogenous PrP dilution linearity (plate 29) ####
  par(mar=c(3,4,3,1))
  
  wt_samples = c(1,2,7,8)
  opt_6d11 %>%
    filter(plate==29) %>%
    mutate(wt_sample = case_when(sample_id %in% wt_samples ~ sample_id,
                                 !(sample_id %in% wt_samples) ~ lag(sample_id))) %>%
    group_by(wt_sample) %>%
    mutate(rel = conc_ngml / conc_ngml[proportion_wt==1]) -> dillin
  
  xlims = ylims = c(0, 1.05)
  xats = seq(0,1,0.05)
  xbigs = unique(dillin$proportion_wt)
  xbiglabs = percent(xbigs)
  yats = seq(0,1,0.05)
  ybigs = unique(dillin$proportion_wt)
  ybiglabs = percent(ybigs)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbiglabs, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='proportion WT brain', cex=0.8)
    axis(side=2, at=yats, labels=NA, tck=-0.025)
    axis(side=2, at=ybigs, labels=NA, tck=-0.05)
    axis(side=2, at=ybigs, labels=ybiglabs, lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='relative PrP concentration', cex=0.8)
  wt_colors = c('1'='#CC0000', '2'='#0099CC', '7'='#33AA44', '8'='#FF9900')
  for (samp in unique(dillin$wt_sample)) {
    dillin %>%
      filter(wt_sample==samp) -> subs
    col = wt_colors[as.character(samp)]
    points(x=subs$proportion_wt, y=subs$rel, pch=19, col=col)
    points(x=subs$proportion_wt, y=subs$rel, type='l', lwd=0.5, col=col)
  }
  abline(a=0, b=1, lty=3)
  mtext('Linear response', side=3, line=1, cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj=0.05, line=0.5); panel = panel + 1

  write_supp_table(dillin %>% select(sample_id, wt_sample, proportion_wt, conc_ngml, rel), 'Endogenous PrP dilution linearity per-sample data, plate 29 (Figure S1A).')

  ### B. ELISA dilution series (plate 30-31) ####
  opt_6d11 %>%
    filter(plate %in% 30:31) -> parallelism

  prop_colors  = c('1'='#003366', '0.15'='#0099CC', '0.05'='#91b2c3', '0'='#c1c1c1')
  dil_offsets  = c('200'=-1, '400'=0, '800'=1)
  group_centers = c('1'=2, '0.15'=6, '0.05'=10, '0'=14)

  parallelism = parallelism %>%
    mutate(xpos  = group_centers[as.character(proportion_wt)] + dil_offsets[as.character(dilution)],
           color = prop_colors[as.character(proportion_wt)])

  parallelism %>%
    group_by(proportion_wt, dilution) %>%
    summarize(.groups='keep',
              n=n(),
              xpos       = first(xpos),
              color      = first(color),
              mean_conc  = mean(conc_ngml),
              l95        = lower(conc_ngml),
              u95        = upper(conc_ngml)) %>%
    ungroup() -> parallelism_smry

  yats  = seq(0, 800, 200)
  xlims = c(0, 16)
  ylims = c(0, 800)
  par(mar=c(3,4,3,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=yats, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=yats, labels=formatC(yats, format='d', big.mark=','), lwd=0, line=-0.5, las=2, cex.axis=0.8)
  mtext(side=2, line=3, text='PrP (ng/mL)', cex=0.7)
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, at=parallelism_smry$xpos, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=1, at=parallelism_smry$xpos, labels=paste0('1:',parallelism_smry$dilution), lwd=0, line=-0.5, cex.axis=0.6, las=2)
  mtext(side=1, at=group_centers, text=c('100%', '15%', '5%', '0%'), line=2, cex=0.8)
  abline(v=c(4, 8, 12), lty=3, col='#AAAAAA')
  suppressWarnings(arrows(x0=parallelism_smry$xpos, y0=parallelism_smry$l95, y1=parallelism_smry$u95,
         code=3, angle=90, length=0.03, lwd=1.5, col=parallelism_smry$color))
  points(x=parallelism_smry$xpos, y=parallelism_smry$mean_conc, pch=19, col=parallelism_smry$color, cex=1.2)
  points(x=parallelism$xpos, y=parallelism$conc_ngml, pch=21, col=parallelism$color, bg='#FFFFFF', cex=0.8)
  mtext('Parallelism', side=3, line=1, cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj=0.05, line=0.5); panel = panel + 1

  write_supp_table(parallelism_smry %>% select(-color), 'Summary of ELISA parallelism dilution series, plates 30-31 (Figure S1B).')

  silence_is_golden = dev.off()
}
### end Figure S1 ####

## Figure 1.  #### 
tell_user('done.\nCreating Figure 1...')
if ('figure-1'=='figure-1') {
  
  resx=300
  png('display_items/figure-1.png',width=6.5*resx,height=4.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,1,1,1,2,2,
                           3,3,3,4,4,4), nrow=2, byrow=T)
  layout(layout_matrix, heights=c(1,1.5))
  
  panel = 1
  
  par(mar=c(3,1,3,2))
  ### A. placeholder / blank for chemical structure or diagram ####
  raster_panel = image_convert(image_read('data/sirna_a_diagram_alny_style.png'),'png')
  plot(as.raster(raster_panel))
  mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.05, line = 0.5); panel = panel + 1
  
  par(mar=c(3,4,3,1))
  ### B. dose-response RNA target engagement ####
  rna_te_meta = tribble(
    ~dose, ~color,
    0,   '#c1c1c1',
    50,  '#91b2c3',
    150, '#0099cc',
    300, '#003366'
  )
  qpcr %>%
    inner_join(rna_te_meta, by='dose') -> rna
  elisa %>%
    filter(plate==33) %>%
    mutate(dose = as.integer(cohort),
           animal = as.integer(sample)) %>%
    mutate(rel = ngml_av / mean(ngml_av[dose==0])) %>%
    inner_join(rna, by=c('dose','animal','expt')) -> rna_protein
  rna_protein %>%
    group_by(dose, color) %>%
    summarize(.groups='keep',
              n=n(),
              rna_mean=mean(rna),
              rna_l95=lower(rna),
              rna_u95=upper(rna),
              prp_mean=mean(rel),
              prp_l95=lower(rel),
              prp_u95=upper(rel)) %>%
    ungroup() -> rna_protein_smry
  xlims = c(0, 1.25)
  ylims = c(0, 1.25)
  barwidth=15
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=1, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=1, at=0:5/4, labels=paste0(formatC(100*0:5/4, format='f',digits=0),'%'), lwd=0, line=-0.5, cex.axis=0.8)
  mtext(side=1, line=1.6, text=expression(residual~italic(Prnp)), cex=0.7)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=paste0(formatC(100*0:5/4, format='f',digits=0),'%'), lwd=0, line=-0.5, las=2, cex.axis=0.8)
  mtext(side=2, line=2.4, text='residual PrP', cex=0.7)
  abline(h=1, lty=3)
  abline(v=1, lty=3)
  barwidth=15
  arrows(x0=rna_protein_smry$rna_l95, x1=rna_protein_smry$rna_u95, y0=rna_protein_smry$prp_mean, code=3, angle=90, length=0.02, lwd=1.5, col=rna_protein_smry$color)
  arrows(x0=rna_protein_smry$rna_mean, y0=rna_protein_smry$prp_l95, y1=rna_protein_smry$prp_u95, code=3, angle=90, length=0.02, lwd=1.5, col=rna_protein_smry$color)
  text(x=rna_protein_smry$rna_mean-0.1, y=rna_protein_smry$prp_mean+0.18, labels=paste0(rna_protein_smry$dose,' µg'), col=rna_protein_smry$color, cex=0.7)
  points(rna_protein$rna, rna_protein$rel, pch=21, col=rna_protein$color, bg='#FFFFFF')
  mtext('ELISA vs. qPCR', side=3, line=1, cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.05, line = 0.5); panel = panel + 1

  write_supp_table(rna_protein_smry %>% select(-color), 'Summary of RNA target engagement dose-response, ELISA vs. qPCR (Figure 1B).')

  ### C. repeat dosing regimens ####
  cohort_meta = tribble(
    ~cohort, ~y, ~color, ~doses,
    'no dose', 6, '#A7A7A7', c(),
    '0 -> 60', 5, '#41ae76', c(0),
    '0, 30 -> 60', 4,'#005824', c(0, 30),
    '0 -> 150', 3, '#c994c7', c(0),
    '0, 30, 120 -> 150', 2, '#91003f', c(0, 30, 120),
  )
  
  cohort_events = tribble(
    ~cohort, ~event, ~day, ~pch,
    'no dose', 'takedown', 60, 15,
    'no dose', 'takedown', 150, 15,
    '0 -> 60', 'dose', 0, 25,
    '0 -> 60', 'takedown', 60, 15,
    '0, 30 -> 60', 'dose', 0, 25,
    '0, 30 -> 60', 'dose', 30, 25,
    '0, 30 -> 60', 'takedown', 60, 15,
    '0 -> 150', 'dose', 0, 25,
    '0 -> 150', 'takedown', 150, 15,
    '0, 30, 120 -> 150', 'dose', 0, 25,
    '0, 30, 120 -> 150', 'dose', 30, 25,
    '0, 30, 120 -> 150', 'dose', 120, 25,
    '0, 30, 120 -> 150', 'takedown', 150, 15
  ) %>%
    inner_join(cohort_meta, by='cohort')
  
  
  cohort_ranges = tribble(
    ~cohort, ~start, ~stop,
    'no dose', 0, 60,
    'no dose', 0, 150,
    '0 -> 60', 0, 60,
    '0, 30 -> 60', 0, 60,
    '0 -> 150', 0, 150,
    '0, 30, 120 -> 150', 0, 150
  ) %>%
    inner_join(cohort_meta, by='cohort')
  
  
  cohort_doses = tribble(
    ~cohort, ~event, ~day, ~pch,
    '0 -> 60', 'dose', 0, 25,
    '0, 30 -> 60', 'dose', 0, 25,
    '0, 30 -> 60', 'dose', 30, 25,
    '0 -> 150', 'dose', 0, 25,
    '0, 30, 120 -> 150', 'dose', 0, 25,
    '0, 30, 120 -> 150', 'dose', 30, 25,
    '0, 30, 120 -> 150', 'dose', 120, 25
  ) %>%
    inner_join(cohort_meta, by='cohort')
  
  
  xlims = c(-30, 160)
  ylims = c(range(cohort_meta$y) + c(-0.5, 0.5))
  par(mar=c(3,8,3,0.5))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  xats = seq(0,150,30)
  axis(side=1, at=xats, labels=NA, cex.axis=.8)
  axis(side=1, at=xats, labels=xats, line=-0.5, lwd=0, cex.axis=.8)
  mtext(side=3, line=1, text='dosing and takedown regimen', cex=0.8)
  # axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  mtext(side=2, at=cohort_meta$y, las=2, text=gsub('->','→',cohort_meta$cohort), cex=0.8)
  points(cohort_doses$day, cohort_doses$y+0.1, pch=cohort_doses$pch, col=cohort_doses$color, bg=cohort_doses$color)
  text(cohort_doses$day, cohort_doses$y+0.2, pos=3, labels='300\nug', col=cohort_doses$color, cex=0.6)
  arrows(x0=cohort_ranges$start, x1=cohort_ranges$stop, y0=cohort_ranges$y, col=cohort_ranges$color, code=2, angle=90, length=0.05)
  mtext(side=1, line=1.6, text='day', cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1
  
  ### D. 6D11 ELISA knockdown ####
  par(mar=c(3,0.5,3,1))
  xlims = c(0, 1.5)
  ylims = c(range(cohort_meta$y) + c(-0.5, 0.5))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:6/4, labels=NA, cex.axis=.8)
  axis(side=1, at=0:6/4, labels=percent(0:6/4), line=-0.5, lwd=0, cex.axis=.8)
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  
  elisa %>%
    filter(expt=='CMR-3191') %>%
    mutate(rel = ngml_av / mean(ngml_av[!is.na(cohort) & cohort=='no dose'], na.rm=T)) %>%
    left_join(cohort_meta, by='cohort') -> elisa_anno
  
  elisa_anno %>% 
    mutate(tenpercent = ngml_av[sample=='MoPosLoQC']/ngml_av[sample=='MoPosHiQC']) %>% 
    select(tenpercent) %>% 
    slice(1) %>% 
    pull() -> tenpercent
  llq = 0.05 * 800 / mean(elisa_anno$ngml_av[elisa_anno$cohort=='no dose'], na.rm=T)
  abline(v=1, lty=3)
  abline(v=llq, lty=2, lwd=0.25)
  mtext(side=3,at=llq,text='LLQ', cex=0.45)
  abline(v=tenpercent, lty=2, lwd=0.25)
  mtext(side=3,at=tenpercent,text='10%\nQC', cex=0.45)
  
  elisa_anno %>%
    group_by(cohort) %>%
    summarize(.groups='keep', 
              n=n(),
              mean = mean(rel), 
              l95 = lower(rel),
              u95 = upper(rel),
              maxr = max(rel)) %>%
    ungroup() %>%
    inner_join(cohort_meta, by='cohort') -> elisa_smry
  barwidth = 0.8
  rect(xleft=rep(0, nrow(elisa_smry)), xright=elisa_smry$mean, ybottom=elisa_smry$y-barwidth/2, ytop=elisa_smry$y+barwidth/2, col=alpha(elisa_smry$color, ci_alpha), border=NA)
  arrows(x0=elisa_smry$l95, x1=elisa_smry$u95, y0=elisa_smry$y, code=3, angle=90, length=0.1, col='#000000')
  points(elisa_anno$rel, elisa_anno$y, pch=21, col=elisa_anno$color, bg='#FFFFFF')
  text(x=elisa_smry$maxr, y=elisa_smry$y, pos=4, labels=percent(elisa_smry$mean, digits=1), cex=0.8)
  mtext(side=3, line=1, text='ELISA', cex=0.8)
  mtext(side=1, line=1.6, text='% residual PrP', cex=0.8)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.05, line = 0.5); panel = panel + 1

  write_supp_table(elisa_smry %>% select(-color,-y,-doses,-maxr), 'Summary of 6D11 ELISA PrP knockdown by cohort (Figure 1D).')

  silence_is_golden = dev.off()
}






## Figure 2.  #### 
tell_user('done.\nCreating Figure 2...')
if ('figure-2'=='figure-2') {
  
  resx=300
  png('display_items/figure-2.png',width=6.5*resx,height=8*resx,res=resx)
  
  layout_matrix = matrix(c(1,1,1,
                           2,3,4,
                           5,6,7,
                           8,8,8), nrow=4, byrow=T)
  layout(layout_matrix)
  
  panel = 1
  
  par(mar=c(3,4,3,10))
  
  a_60_dpi_surgeries = c(60,90,179,266,361,447)
  
  ### A. 60 dpi survival ####
  sf = make_sf(survdata %>% filter(expt=='A' & grepl('(60|uninoc)',cohort)))
  sv = make_sv(survdata %>% filter(expt=='A' & grepl('(60|uninoc)',cohort) & grepl('(aCSF|sirnaA)',cohort)))
  sv_pval(sv)
  plot_sf(sf, anno_dpi=a_60_dpi_surgeries)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.09, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='survival - RML 60 dpi pre-symptomatic treatment', cex=0.8)

  write_supp_table(sf_table(sf), 'Survival summary for RML 60 dpi pre-symptomatic treatment (Figure 2A).')

  ### B. 60 dpi weights ####
  par(mar=c(3,4,3,1))
  baseline_dpi = 73
  nw = norm_wts(weights_raw=weights,
                    baseline_dpi=73,
                    main_tbl = survdata %>% filter(expt=='A' & grepl('(60|uninoc)',cohort)),
                    meta_tbl = meta %>% filter(!exclude))
  plot_wts(nw)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='weights', cex=0.8)

  write_supp_table(nw[[2]] %>% ungroup() %>% select(-color,-lty), 'Summary of weight change for RML 60 dpi pre-symptomatic treatment (Figure 2B).')
  
  ### C. 60 dpi NfL ####
  
  nfl %>%
    inner_join(blind, by='animal') %>%
    inner_join(meta %>% filter(!exclude), by='cohort') %>%
    filter(grepl('60',cohort) | grepl('uninoc',cohort)) -> nfl_anno
  nfl_anno %>%
    group_by(dpi, cohort, color, lty) %>%
    summarize(.groups='keep', 
              n=n(),
              mean=mean(conc), 
              l95=lower(conc),
              u95=upper(conc)) %>%
    ungroup() -> nfl_smry
  
  nfl_smry$ks_p = as.numeric(NA)
  for (i in 1:nrow(nfl_smry)) {
    dpi = nfl_smry$dpi[i]
    this_rows = nfl_anno$dpi == dpi & nfl_anno$cohort == nfl_smry$cohort[i]
    ctrl_rows = nfl_anno$dpi == dpi & nfl_anno$cohort == 'A-uninoc'
    ks_obj = suppressWarnings(ks.test(nfl_anno$conc[this_rows], nfl_anno$conc[ctrl_rows], alternative='two.sided'))
    nfl_smry$ks_p[i] = ks_obj$p.value
  }
  nfl_smry$p_symb = ''
  nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.05] = '*'
  nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.01] = '**'
  nfl_smry$p_symb[!is.na(nfl_smry$ks_p) & nfl_smry$ks_p < 0.001] = '***'
  
  par(mar=c(3,4,3,1))
  xlims=c(0,500)
  xats = seq(0,500,50)
  xbigs = seq(0,500,100)
  ylims = c(0, 1500)
  yats = seq(0, 1500, 500)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xbigs, labels=NA, tck=-0.05)
  axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
  axis(side=1, at=xats, labels=NA, tck=-0.025)
  mtext(side=1, line=1.5, text='days post-inoculation', cex=0.8)
  axis(side=2, at=yats, labels=NA, tck=-0.05)
  axis(side=2, at=yats, labels=formatC(yats,big.mark=','), lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='NfL (pg/mL)', cex=0.8)
  par(xpd=T)
  for (coh in unique(nfl_smry$cohort)) {
    subs = subset(nfl_smry, cohort==coh)
    color = subs$color[1]
    points(subs$dpi, subs$mean, type='l', lwd=3, col=color)
    points(subs$dpi, subs$mean, pch=15, col=color)
    polysubs = subset(subs, !is.na(l95)) 
    polygon(x=c(polysubs$dpi,rev(polysubs$dpi)), y=c(polysubs$u95,rev(polysubs$l95)), col=alpha(color,ci_alpha), border=NA)
  }
  points(nfl_anno$dpi, nfl_anno$conc, pch=21, col=nfl_anno$color, bg='#FFFFFF')
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='plasma NfL', cex=0.8)

  write_supp_table(nfl_smry %>% select(-color,-lty), 'Summary of plasma NfL for RML 60 dpi pre-symptomatic treatment (Figure 2C).')

  ### D. 60 dpi nests ####
  nn = norm_nests(nests,
                main_tbl = survdata %>% filter(expt=='A' & grepl('(60|uninoc)',cohort)))
  par(mar=c(3,4,3,1))
  plot_nests(nn)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='nests', cex=0.8)

  nn %>%
    group_by(cohort, disp, nestdpi) %>%
    summarize(.groups='keep', n=n(), l95=lower(mean), u95=upper(mean), mean=mean(mean)) %>%
    ungroup() -> nest_smry_60dpi
  write_supp_table(nest_smry_60dpi, 'Summary of nest scores for RML 60 dpi pre-symptomatic treatment (Figure 2D).')

  ### E. 140 dpi survival ####
  sf = make_sf(survdata %>% filter(expt=='A' & grepl('(140)',cohort)))
  sv = make_sv(survdata %>% filter(expt=='A' & grepl('(140|uninoc)',cohort) & grepl('(aCSF|sirnaA)',cohort)))
  sv_pval(sv)
  plot_sf(sf, xlims=c(0,250), anno_dpi=c(140,174), show_legend = F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=0.6, text='\nsurvival - RML 140 dpi\nlate treatment', cex=0.8)

  write_supp_table(sf_table(sf), 'Survival summary for RML 140 dpi late treatment (Figure 2E).')

  nw = norm_wts(weights_raw=weights,
                baseline_dpi=73, 
                main_tbl = survdata %>% filter(expt=='A' & grepl('(140|uninoc)',cohort)),
                meta_tbl = meta %>% filter(!exclude))
  
  ### F. 126 dpi weight evidence of symptoms ####
  
  nw = norm_wts(weights_raw=weights, 
                baseline_dpi=c(75,77), 
                main_tbl = survdata %>% filter(expt=='YA'),
                meta_tbl = meta %>% filter(!exclude))
  
  
  uninoc_col = '#A9A9A9'
  rml_col = '#8E2323'
  meta_125 = tibble(x=1:2, color =  c(uninoc_col, rml_col), inoculum=c('none','RML'))
  
  
  nw[[1]] %>%
    filter(wtdpi == 125) %>%
    inner_join(blind %>% filter(expt=='YA') %>% distinct(cohort, inoculum), by='cohort') %>%
    select(inoculum, rel) %>%
    inner_join(meta_125, by='inoculum') -> wt_125
  wt_125 %>%
    group_by(inoculum, x, color) %>%
    summarize(.groups='keep',
              n=n(),
              mean=mean(rel),
              l95=lower(rel),
              u95=upper(rel)) %>%
    ungroup() -> wt_125_smry
  
  
  xlims=c(0.5, 2.5)
  ylims = c(-0.25,0.4)
  yats = seq(-0.3,0.5,0.1)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
  mtext(side=1, line=0.25, at=meta_125$x, text=meta_125$inoculum, cex=0.8)
  axis(side=2, at=yats, labels=NA, tck=-0.05)
  axis(side=2, at=yats, labels=percent(yats), lwd=0, line=-0.5, las=2)
  mtext(side=2, line=2.75, text='weight (% change)', cex=0.8)
  abline(h=0)
  weights_uninoc = wt_125$rel[wt_125$inoculum=='none']
  weights_rml = wt_125$rel[wt_125$inoculum=='RML']
  weights_t_obj = t.test(weights_uninoc,weights_rml)
  set.seed(1)
  points(x=jitter(rep(1,length(weights_uninoc)),amount=.25), y=weights_uninoc, pch=21, col=uninoc_col, bg='#FFFFFF')
  points(x=jitter(rep(2,length(weights_rml)),amount=.25), y=weights_rml, pch=21, col=rml_col, bg='#FFFFFF')
  barwidth=0.3
  segments(x0=wt_125_smry$x-barwidth, x1=wt_125_smry$x+barwidth, y0=wt_125_smry$mean, col=wt_125_smry$color, lwd=2)
  arrows(x0=wt_125_smry$x, y0=wt_125_smry$l95, y1=wt_125_smry$u95, code=3, angle=90, length=0.05, col=wt_125_smry$color, lwd=2)
   mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='weights at 125 dpi', cex=0.8)

  write_supp_table(wt_125_smry %>% select(-color), 'Summary of weight change at 125 dpi, uninoculated vs. RML (Figure 2F).')

  ### G. 126 dpi nest evidence of symptoms ####
  nn = norm_nests(nests,
                  main_tbl = survdata %>% filter(expt=='YA'))
  
  nn %>%
    filter(nestdpi == 125) %>%
    inner_join(blind %>% filter(expt=='YA') %>% distinct(cohort, inoculum), by='cohort') %>%
    select(cage, inoculum, mean) %>%
    inner_join(meta_125, by='inoculum') -> nest_125
  nest_125 %>%
    rename(nest_score=mean) %>%
    group_by(inoculum, x, color) %>%
    summarize(.groups='keep',
              n=n(),
              mean=mean(nest_score),
              l95=lower(nest_score),
              u95=upper(nest_score)) %>%
    ungroup() ->  nest_125_smry
  xlims = c(0.5,2.5)
  ylims = c(0,2.05)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
  mtext(side=1, line=0.25, at=meta_125$x, text=meta_125$inoculum, cex=0.8)
  axis(side=2, at=0:2, las=2)
  axis(side=2, at=0:4/2, las=2, labels=NA, tck=-0.02)
  mtext(side=2, line=2.75, text='nest score', cex=0.8)
  nests_uninoc = nest_125$mean[nest_125$inoculum=='none']
  nests_rml = nest_125$mean[nest_125$inoculum=='RML']
  nests_t_obj = t.test(nests_uninoc,nests_rml)
  set.seed(1)
  points(x=jitter(rep(1,length(nests_uninoc)),amount=.25), y=nests_uninoc, pch=21, col=uninoc_col, bg='#FFFFFF')
  points(x=jitter(rep(2,length(nests_rml)),amount=.25), y=nests_rml, pch=21, col=rml_col, bg='#FFFFFF')
  barwidth=0.3
  segments(x0=nest_125_smry$x-barwidth, x1=nest_125_smry$x+barwidth, y0=nest_125_smry$mean, col=nest_125_smry$color, lwd=2)
  suppressWarnings(arrows(x0=nest_125_smry$x, y0=nest_125_smry$l95, y1=nest_125_smry$u95, code=3, angle=90, length=0.05, col=nest_125_smry$color, lwd=2))
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='nests at 125 dpi', cex=0.8)

  write_supp_table(nest_125_smry %>% select(-color), 'Summary of nest scores at 125 dpi, uninoculated vs. RML (Figure 2G).')

  dpi_125_t_tests = tribble(
    ~outcome, ~t_test_pval, ~mean_rml, ~mean_uninoc, ~n_rml, ~n_uninoc,
    'weights', weights_t_obj$p.value, weights_t_obj$estimate[1], weights_t_obj$estimate[2], length(weights_rml), length(weights_uninoc),
    'nests', nests_t_obj$p.value, nests_t_obj$estimate[1], nests_t_obj$estimate[2], length(nests_rml), length(nests_uninoc),
  )
  write_supp_table(dpi_125_t_tests, 'T tests for 125 dpi weights and nests.')
  
  ### H. 126 dpi survival ####
  par(mar=c(3,4,3,10))
  ya_126_surgeries = c(126,156,247,338,429)
  sf = make_sf(survdata %>% filter(expt=='YA'))
  plot_sf(sf, anno_dpi=ya_126_surgeries, show_legend=F)
  leg = meta %>% filter(grepl('^YA',cohort) & !exclude)
  disp_order = c('aCSF','late 50 ug once','late 150 ug once','late 300 ug once','late chronic 300 ug','uninoculated')
  leg = leg[match(disp_order, leg$disp),]
  par(xpd=T)
  legend(x=500, y=1, bty='n', cex=.8, lwd=2, legend=leg$disp,
         col=leg$color, lty=leg$lty,
         text.col=leg$color)
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.09, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='survival - RML 126 dpi symptomatic treatment', cex=0.8)

  write_supp_table(sf_table(sf), 'Survival summary for RML 126 dpi symptomatic treatment (Figure 2H).')

  ya_overall_sdiff   = make_sv(survdata %>% filter(expt=='YA'))
  sv_pval(ya_overall_sdiff)
  ya_050_once_sdiff  = make_sv(survdata %>% filter(expt=='YA', disp %in% c('aCSF', 'late 50 ug once')))
  ya_150_once_sdiff  = make_sv(survdata %>% filter(expt=='YA', disp %in% c('aCSF', 'late 150 ug once')))
  ya_300_once_sdiff  = make_sv(survdata %>% filter(expt=='YA', disp %in% c('aCSF', 'late 300 ug once')))
  ya_300_chron_sdiff = make_sv(survdata %>% filter(expt=='YA', disp %in% c('aCSF', 'late chronic 300 ug')))
  silence_is_golden = dev.off()
}


## Figure S2.  ####
tell_user('done.\nCreating Figure S2...')
if ('figure-s2'=='figure-s2') {

  resx=300
  png('display_items/figure-s2.png',width=6.5*resx,height=4.5*resx,res=resx)

  layout_matrix = matrix(c(1,2), nrow=1, byrow=T)
  layout(layout_matrix)

  panel = 1

  ### A. YA weights ####
  par(mar=c(3,4,3,1))
  nw = norm_wts(weights_raw=weights,
                baseline_dpi=c(75,77),
                main_tbl = survdata %>% filter(expt=='YA'),
                meta_tbl = meta %>% filter(!exclude))
  ya_wt_ylims = c(-0.3, 0.6)
  ya_wt_yats = seq(-0.3, 0.6, 0.1)
  plot_wts(nw, anno_dpi=ya_126_surgeries, ylims=ya_wt_ylims, yats=ya_wt_yats, show_legend=F)
  leg = meta %>% filter(grepl('^YA',cohort) & !exclude)
  disp_order = c('aCSF','late 50 ug once','late 150 ug once','late 300 ug once','late chronic 300 ug','uninoculated')
  leg = leg[match(disp_order, leg$disp),]
  legend(x=0, y=0.5, bty='n', cex=.7, lwd=2, legend=leg$disp,
         col=leg$color, lty=leg$lty,
         text.col=leg$color)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='weights', cex=0.8)

  write_supp_table(nw[[2]] %>% ungroup() %>% select(-color,-lty), 'Summary of weight change for RML 126 dpi symptomatic treatment (Figure S2A).')

  ### B. YA nests ####
  nn = norm_nests(nests,
                  main_tbl = survdata %>% filter(expt=='YA'))
  par(mar=c(3,4,3,1))
  plot_nests(nn, anno_dpi=ya_126_surgeries, show_legend=F)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='nests', cex=0.8)

  nn %>% 
    mutate(approx_dpi = case_when(cage %in% c('YAJ','YAU')  ~ nestdpi + 5,
                                  TRUE ~ nestdpi)) %>%
    # above code handles the fact that inoc dates are staggered 5 days so 
    # they don't group properly and you get zigzag effects. note this issue
    # does not affect the plot because that is gaussian smoothed anyway
    group_by(cohort, disp, approx_dpi) %>%
    summarize(.groups='keep', n=n(), mean=mean(mean), l95=lower(mean), u95=upper(mean)) %>%
    ungroup() -> nest_smry_ya
  write_supp_table(nest_smry_ya, 'Summary of nest scores for RML 126 dpi symptomatic treatment (Figure S2B).')

  silence_is_golden = dev.off()
}
### end Figure S2 ####



## Figure 4.  #### 
tell_user('done.\nCreating Figure 4...')
if ('figure-4'=='figure-4') {
  
  resx=300
  png('display_items/figure-4.png',width=6.5*resx,height=4.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,1,
                           2,3), nrow=2, byrow=T)
  layout(layout_matrix)
  
  panel = 1
  
  par(mar=c(3,4,3,10))

  
  ### A. AF survival ####
  af_surgeries = c(56, 86, 189, 266, 357, 448)
  sf = make_sf(survdata %>% filter(expt=='AF'))
  sv_het = make_sv(survdata %>% filter(expt=='AF' & grepl('(HET|WT-NoTx)',cohort)))
  sv_pval(sv_het)
  sv_si = make_sv(survdata %>% filter(expt=='AF' & grepl('(WT-siRNAA|WT-NoTx)',cohort)))
  sv_pval(sv_si)
  
  plot_sf(sf, anno_dpi=af_surgeries)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.09, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='survival - Fukuoka-1 56 dpi pre-symptomatic treatment', cex=0.8)

  write_supp_table(sf_table(sf), 'Survival summary for Fukuoka-1 56 dpi pre-symptomatic treatment (Figure 4A).')

  ### B. AF weights ####
  par(mar=c(3,4,3,1))
  nw = norm_wts(weights_raw=weights,
                baseline_dpi=66,
                main_tbl = survdata %>% filter(expt=='AF'),
                meta_tbl = meta %>% filter(!exclude))
  plot_wts(nw)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='weights', cex=0.8)

  write_supp_table(nw[[2]] %>% ungroup() %>% select(-color,-lty), 'Summary of weight change for Fukuoka-1 56 dpi pre-symptomatic treatment (Figure 4B).')

  ### C. AF nests ####
  nn = norm_nests(nests,
                  main_tbl = survdata %>% filter(expt=='AF'))
  par(mar=c(3,4,3,1))
  plot_nests(nn)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='nests', cex=0.8)

  nn %>%
    group_by(cohort, disp, nestdpi) %>%
    summarize(.groups='keep', n=n(), l95=lower(mean), u95=upper(mean), mean=mean(mean)) %>%
    ungroup() -> nest_smry_af
  write_supp_table(nest_smry_af, 'Summary of nest scores for Fukuoka-1 56 dpi pre-symptomatic treatment (Figure 4C).')


  silence_is_golden = dev.off()
}


## Figure 5.  #### 
tell_user('done.\nCreating Figure 5...')
if ('figure-5'=='figure-5') {
  
  resx=300
  png('display_items/figure-5.png',width=6.5*resx,height=4.5*resx,res=resx)
  
  layout_matrix = matrix(c(1,1,1,1,2,2,2,2,
                           3,3,3,3,3,3,4,4), nrow=2, byrow=T)
  layout(layout_matrix, heights=c(1,1.25))
  
  panel = 1
  
  par(mar=c(7,4,3,1))
  
  rtquic_titers %>%
    rename(animal = sample) %>%
    inner_join(main, by='animal') %>%
    inner_join(blind %>% select(-cage), by=c('animal','expt')) %>%
    inner_join(meta, by='cohort') %>%
    select(expt, animal, titer, dpi, inoculum, treatment, cohort, color, disp) -> rtq_titer_anno
  
  
  plot_rtquic = function(rtq, 
                         xlims=c(0,500), 
                         #xats = seq(0,500,50),
                         #xbigs = seq(0,500,100),
                         ylims = c(5, 13),
                         yats = log10(rep(1:9, 8) * 10^rep(5:12, each=9)),
                         ybigs = 5:12,
                         llq = 3 + -log10(0.002),
                         ci_mode = 'polygon', # other option is arrows
                         show_legend=F) {
    
    plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
    #axis(side=1, at=xbigs, labels=NA, tck=-0.05)
    #axis(side=1, at=xbigs, lwd=0, labels=xbigs, line=-0.5)
    #axis(side=1, at=xats, labels=NA, tck=-0.025)
    
    axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
    axis(side=2, at=yats, tck=-0.025, labels=NA)
    axis(side=2, at=ybigs, tck=-0.05, labels=ybigs, las=2, cex.axis=0.75)
    
    
    rtq %>%
      mutate(grouping = paste0(inoculum, '-', treatment)) %>%
      group_by(grouping, inoculum, treatment, x, color) %>%
      summarize(.groups='keep',
                n=n(),
                mean = mean(titer),
                l95 = lower(titer),
                u95 = upper(titer)) %>%
      ungroup() -> rtq_smry
    
    barwidth = diff(xlims)/50
    
    if (ci_mode=='polygon') {
      for (this_group in unique(rtq_smry$grouping)) {
        subs = rtq_smry %>% filter(grouping==this_group)
        points(x=subs$x, y=subs$mean, col=subs$color, lwd=1.5, type='l')
        polygon(x=c(subs$x,rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
      }
    } else if (ci_mode == 'arrows') {
      segments(x0=rtq_smry$x-barwidth, x1=rtq_smry$x+barwidth, y0=rtq_smry$mean, col=rtq_smry$color, lwd=2)
      suppressWarnings(arrows(x0=rtq_smry$x, y0=rtq_smry$l95, y1=rtq_smry$u95, code=3, angle=90, length=0.05, col=rtq_smry$color, lwd=2))
    }
    
    set.seed(1)
    rtq$xdisp = jitter(rtq$x, amount=barwidth)
    points(rtq$xdisp, rtq$titer, pch=21, col=rtq$color, bg='#FFFFFF')
    
    abline(h=llq, lty=3)
    mtext(side=4, line=0.25, at=llq, text='LLQ', cex=0.4, las=2)
    
  }
  
  
  
  tranche_plotter = function(xbars, baseline=0, tck=0.01, lwd=1, colorcol=NA) {
    for (i in 1:nrow(xbars)) {
      if (is.na(colorcol)) {
        color='#000000'
      } else {
        color = xbars[i,colorcol] %>% pull()
      }
      axis(side=1, line=baseline, at=c(xbars$xmin[i],xbars$xmax[i]), tck=tck, labels=NA, lwd=lwd, col=color)
      mtext(side=1, line=baseline-0.15, padj=1, at=xbars$x[i], text=xbars$xlabels[i], cex=0.45)
    }
  }
  
  
  ### A. QS RML ####
  par(mar=c(3,4,3,1))
  
  rtq_titer_anno %>%
    filter(expt=='QS') %>%
    filter(inoculum=='RML') %>%
    mutate(x=dpi) -> rtq
  
  plot_rtquic(rtq, xlims=c(0,250))
  mtext(side=2, line=2.25, text='log10 seeds/g', cex=0.7)
  barwidth=8
  xbars = tibble(
    x       = c(60, 90, 120, 150, 240),
    xlabels = paste0(c(60, 90, 120, 150, 240),'\ndpi'),
  ) %>% mutate(xmin = x - barwidth, xmax = x + barwidth)
  tranche_plotter(xbars)
  
  rtq_titer_anno %>%
    filter(expt=='QS') %>%
    distinct(treatment,color) %>%
    mutate(treatment = gsub('het','Het KO', treatment)) %>%
    arrange(desc(treatment)) -> leg
  
  legend('topright', leg$treatment, col=leg$color, pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, lwd=1.5, bty='n')
  
  
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='RML natural history', cex=0.7)
  
  

  ### C. QS Fukuoka-1 ####
  
  rtq_titer_anno %>%
    filter(expt=='QS') %>%
    filter(inoculum=='Fukuoka-1') %>%
    mutate(x=dpi)  -> rtq
  
  
  par(mar=c(3,4,3,1))
  plot_rtquic(rtq, xlims=c(0,250))
  mtext(side=2, line=2.25, text='log10 seeds/g', cex=0.7)
  barwidth=8
  xbars = tibble(
    x       = c(60, 90, 120, 150, 240),
    xlabels = paste0(c(60, 90, 120, 150, 240),'\ndpi'),
  ) %>% mutate(xmin = x - barwidth, xmax = x + barwidth)
  tranche_plotter(xbars)
  
  legend('topright', leg$treatment, col=leg$color, pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, lwd=1.5, bty='n')
  
  
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  mtext(side=3, line=1, text='Fukuoka-1 natural history', cex=0.7)
  
  rtq_titer_anno %>%
    filter(expt=='QS') -> rtq_model_data
  
  rtq_model_data %>%
    group_by(dpi, inoculum, treatment) %>%
    summarize(.groups='keep',
              n=n(),
              mean=mean(titer),
              l95 = lower(titer),
              u95 = upper(titer)) %>%
    ungroup() %>%
    rename(genotype=treatment) %>%
    arrange(inoculum, dpi, genotype) -> rtq_nathist_smry
  
  m = lm(titer ~ dpi + inoculum * treatment, data=rtq_model_data)
  summary(m)$coefficients %>%
    as_tibble(rownames='variable') %>%
    clean_names() %>%
    rename(p = pr_t) -> nathist_lm
  
  write_supp_table(nathist_lm, 'Linear model results for natural history of prion titer of RML and Fukuoka-1 in WT and Het KO mice.')
  
  ### B. RML siRNA effect  ####
  
  a_ya_ar_meta = tribble(
    ~cohort, ~x, ~disp_layer2, ~disp_layer1,
    'AR-untreated-60',  2, '60 dpi',  'untreated',
    'AR-untreated-120', 3, '120 dpi', 'untreated',
    'AR-treated-120',   3, '120 dpi', 'siRNA A',
    'AR-treated-210',   4, '210 dpi', 'siRNA A',
    'A-60-aCSF',  5, 'terminal', 'aCSF',
    'A-140-aCSF', 5, 'terminal', 'aCSF',
    'YA-0-once',  5, 'terminal', 'aCSF',
    'A-60-sirnaB',  6, 'terminal', 'siRNA B', 
    'A-140-sirnaB', 6, 'terminal', 'siRNA B', 
    'YA-50-once', 7, 'terminal', 'siRNA A\n50 µg\n1x', 
    'YA-150-once', 8, 'terminal', 'siRNA A\n150 µg\n1x', 
    'YA-300-once', 9, 'terminal', 'siRNA A\n300 µg\n1x', 
    'YA-300-chronic', 10, '500 dpi', 'siRNA A\n300 µg\nlate\nchronic',
    'A-60-sirnaA', 11, '500 dpi', 'siRNA A\n300 µg\nearly\nchronic', 
    'A-uninoc', 12, '500 dpi', 'uninoc-\nulated', 
    'YA-uninoc', 12, '500 dpi', 'uninoc-\nulated',
  ) %>% 
    inner_join(meta %>% select(cohort, color), by='cohort')

  rtq_titer_anno %>%
    select(-color) %>%
    filter(expt %in% c('AR','A','YA')) %>%
    inner_join(a_ya_ar_meta, by='cohort') -> rtq
  
  rtq %>%
    mutate(treatment_group = gsub('\n',' ',disp_layer1),
           timepoint = disp_layer2) %>%
    group_by(x, timepoint, treatment_group, color, inoculum) %>%
    summarize(.groups='keep',
              n=n(),
              mean=mean(titer),
              l95 = lower(titer),
              u95 = upper(titer)) %>%
    ungroup() %>%
    arrange(x) -> rtq_rml_sirna_smry
  
  # AR - 120
  t.test(rtq_titer_anno$titer[rtq_titer_anno$cohort=='AR-untreated-120'], rtq_titer_anno$titer[rtq_titer_anno$cohort=='AR-treated-120'])
  # range of mean terminal titers
  range(rtq_rml_sirna_smry$mean[rtq_rml_sirna_smry$timepoint=='terminal'])
  # range of individual terminal titers
  range(rtq$titer[rtq$disp_layer2=='terminal'])
  anova_among_terminal = aov(titer ~ disp_layer1, data= rtq %>% filter(disp_layer2=='terminal'))
  summary(anova_among_terminal)
  # A & YA terminal vs. late chronic
  t.test(rtq$titer[rtq$disp_layer2=='terminal'], rtq$titer[rtq$cohort=='YA-300-chronic'])
  t.test(rtq$titer[rtq$disp_layer2=='terminal'], rtq$titer[rtq$cohort=='A-60-sirnaA'])
  
  par(mar=c(7,4,3,1))
  plot_rtquic(rtq, xlims=c(1.5,12.5), ci_mode='arrows')
  mtext(side=2, line=2.25, text='log10 seeds/g', cex=0.7)
  abline(v=4.5, lwd=.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  barwidth = 0.25
  a_ya_ar_meta %>%
    filter(x > 4) %>%
    group_by(disp_layer1, color) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer1) -> layer1
  tranche_plotter(layer1, baseline=0, tck=0.00, colorcol='color', lwd=5)
  
  a_ya_ar_meta %>%
    group_by(disp_layer2) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer2) -> layer2
  tranche_plotter(layer2, baseline=4.5, tck=0.01)
  
  a_ya_ar_meta %>%
    filter(x < 4.5 & disp_layer1 == 'untreated') %>%
    group_by(disp_layer1, color) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer1) -> layer1
  tranche_plotter(layer1, baseline=-0.1, tck=0.00, colorcol='color', lwd=5)
  
  
  a_ya_ar_meta %>%
    filter(x < 4.5 & disp_layer1 == 'siRNA A') %>%
    group_by(disp_layer1, color) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer1) -> layer1
  tranche_plotter(layer1, baseline=0.1, tck=0.00, colorcol='color', lwd=5)
  
  mtext(side=1, line=3, at=rtq_rml_sirna_smry$x, text=paste0('N=',rtq_rml_sirna_smry$n), cex=0.4)
  
  mtext(side=3, line=1, at=3, text='RML + siRNA\ntimed sacrifice', cex=0.7)
  mtext(side=3, line=1, at=8.5, text='RML + siRNA\nsurvival studies', cex=0.7)
  
  
  
  ### D. Fukuoka-1 siRNA effect ####
  
  af_s_meta = tribble(
    ~cohort, ~x, ~disp_layer2, ~disp_layer1, ~color,
    'S-Fu1-PBS',  2, 'terminal', 'PBS', '#E04006',
    'AF-WT-siRNAA', 3, '500 dpi', 'siRNA A\n300 µg\nearly\nchronic', '#003366'
  )
  
  rtq_titer_anno %>%
    select(-color) %>%
    filter(expt %in% c('AF','S')) %>%
    inner_join(af_s_meta, by='cohort') -> rtq
  
  
  rtq %>%
    group_by(x, color, inoculum, treatment, disp_layer2) %>%
    summarize(.groups='keep',
              n=n(),
              mean=mean(titer),
              l95 = lower(titer),
              u95 = upper(titer)) %>%
    ungroup() %>%
    rename(timepoint=disp_layer2) %>%
    arrange(x) -> rtq_fu1_sirna_smry
  
  bind_rows(
    rtq_nathist_smry %>% mutate(panel='Fig 5A-B: RML/Fukuoka-1 natural history'),
    rtq_rml_sirna_smry %>% select(-color) %>% mutate(panel='Fig 5C: RML + siRNA'),
    rtq_fu1_sirna_smry %>% select(-color) %>% mutate(panel='Fig 5D: Fukuoka-1 + siRNA')
  ) %>%
    mutate(timepoint = ifelse(!is.na(dpi), paste0(dpi, ' dpi'), timepoint),
           treatment = coalesce(treatment_group, treatment, 'none'),
           genotype = coalesce(genotype, 'WT')) %>%
    select(panel, timepoint, inoculum, genotype, treatment, x, n, mean, l95, u95) -> rtq_all_panels_smry

  write_supp_table(rtq_all_panels_smry, 'Summary of RT-QuIC titer data for RML and Fukuoka-1 natural history and siRNA treatment studies (Figure 5A-D).')
  
  t.test(rtq$titer[rtq$treatment=='PBS'], rtq$titer[rtq$treatment=='siRNA-A'])
  
  par(mar=c(7,4,3,1))
  plot_rtquic(rtq, xlims=c(1.5,3.5), ci_mode='arrows')
  mtext(side=2, line=2.25, text='log10 seeds/g', cex=0.7)
  abline(v=4.5, lwd=.5)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.11, line = 0.5); panel = panel + 1
  
  barwidth = 0.25
  af_s_meta %>%
    group_by(disp_layer1, color) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer1) -> layer1
  tranche_plotter(layer1, baseline=0, tck=0.00, colorcol='color', lwd=5)
  
  af_s_meta %>%
    group_by(disp_layer2) %>%
    summarize(.groups='keep', 
              xmin = min(x)-barwidth,
              xmax = max(x)+barwidth) %>%
    ungroup() %>%
    mutate(x = (xmin+xmax)/2) %>%
    rename(xlabels = disp_layer2) -> layer2
  tranche_plotter(layer2, baseline=4.5, tck=0.01)
  
  mtext(side=1, line=3, at=rtq_fu1_sirna_smry$x, text=paste0('N=',rtq_fu1_sirna_smry$n), cex=0.4)
  
  mtext(side=3, line=1, text='Fukuoka-1 + siRNA\nsurvival study', cex=0.7)

  silence_is_golden = dev.off()
}


## Figure 6.  ####
tell_user('done.\nCreating Figure 6...')
if ('figure-6'=='figure-6') {

  resx=300
  png('display_items/figure-6.png',width=6.5*resx,height=5*resx,res=resx)

  layout_matrix = matrix(c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),
                           rep(6,8),rep(7,11),rep(8,11)), nrow=2, byrow=T)
  layout(layout_matrix, heights=c(1,1.25))

  panel = 1

  mouse_dose_meta = tibble(
    dose_mg   = c(0, 0.0375, 0.075, 0.15, 0.3),
    dose_disp = c('0','37.5','75','150','300'),
    x         = 1:5,
    color     = colorRampPalette(c('#c1c1c1','#0099cc','#003366'))(5)
  )

  nhp_dose_meta = tibble(
    dose_mg   = c(0, 10, 60),
    dose_disp = c('0 mg','10 mg','60 mg'),
    color     = c('#A9A9A9','#33CCFF','#003366'),
    offset    = c(-0.25, 0, 0.25)
  )

  region_meta = tibble(
    tissue      = c('Prefrontal Cortex','Hippocampus','Thalamus','Hypothalamus'),
    tissue_disp = c('prefrontal\ncortex','hippo-\ncampus','thalamus','hypo-\nthalamus'),
    region_x    = 1:4
  )

  mouse_barwidth = 0.35


  ### A. mouse ICV dose response, mRNA, 28 days ####
  par(mar=c(3,4,3,1))

  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='mRNA', days_post_dose==28) %>%
    inner_join(mouse_dose_meta, by='dose_mg') -> panel_a_data

  panel_a_data %>%
    group_by(dose_mg, dose_disp, x, color) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() %>%
    arrange(x) -> panel_a_smry

  ylims = c(0, 1.3)
  xlims_a = c(0.5,5.5)
  plot(NA, NA, xlim=xlims_a, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims_a, labels=NA, lwd.ticks=0)
  mtext(side=1, at=panel_a_smry$x, text=panel_a_smry$dose_disp, las=2, line=0.1, padj=0, cex=0.4)
  mtext(side=1, line=1.6, text='dose (µg), 28d', cex=0.7)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text=expression(residual~italic(PRNP)), cex=0.7)
  abline(h=1, lty=2)
  segments(x0=panel_a_smry$x-mouse_barwidth, x1=panel_a_smry$x+mouse_barwidth, y0=panel_a_smry$mean, col=panel_a_smry$color, lwd=2)
  suppressWarnings(arrows(x0=panel_a_smry$x, y0=panel_a_smry$l95, y1=panel_a_smry$u95, code=3, angle=90, length=0.05, col=panel_a_smry$color, lwd=1.5))
  set.seed(1)
  points(jitter(panel_a_data$x, amount=0.15), panel_a_data$residual, pch=21, col=panel_a_data$color, bg='#FFFFFF')
  mtext(side=3, line=0.5, text='mRNA', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.15, line=0.2); panel = panel + 1


  ### B. mouse ICV dose response, protein, 28 days ####
  par(mar=c(3,4,3,1))

  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='protein', days_post_dose==28) %>%
    inner_join(mouse_dose_meta, by='dose_mg') -> panel_b_data

  panel_b_data %>%
    group_by(dose_mg, dose_disp, x, color) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() %>%
    arrange(x) -> panel_b_smry

  plot(NA, NA, xlim=xlims_a, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims_a, labels=NA, lwd.ticks=0)
  mtext(side=1, at=panel_b_smry$x, text=panel_b_smry$dose_disp, las=2, line=0.1, padj=0, cex=0.4)
  mtext(side=1, line=1.6, text='dose (µg), 28d', cex=0.7)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text='residual PrP', cex=0.7)
  abline(h=1, lty=2)
  segments(x0=panel_b_smry$x-mouse_barwidth, x1=panel_b_smry$x+mouse_barwidth, y0=panel_b_smry$mean, col=panel_b_smry$color, lwd=2)
  suppressWarnings(arrows(x0=panel_b_smry$x, y0=panel_b_smry$l95, y1=panel_b_smry$u95, code=3, angle=90, length=0.05, col=panel_b_smry$color, lwd=1.5))
  set.seed(1)
  points(jitter(panel_b_data$x, amount=0.15), panel_b_data$residual, pch=21, col=panel_b_data$color, bg='#FFFFFF')
  mtext(side=3, line=0.5, text='PrP', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.15, line=0.2); panel = panel + 1


  ### C. mouse mRNA & protein correlation, 28 days ####
  par(mar=c(3,2.5,3,1))

  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='mRNA', days_post_dose==28) %>%
    select(animal_id, dose_mg, mrna=residual) -> mrna_28
  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='protein', days_post_dose==28) %>%
    select(animal_id, prot=residual) -> prot_28

  mrna_28 %>%
    inner_join(prot_28, by='animal_id') %>%
    inner_join(mouse_dose_meta, by='dose_mg') -> panel_c_data

  panel_c_data %>%
    group_by(dose_mg, dose_disp, color) %>%
    summarize(.groups='keep', n=n(), mrna_mean=mean(mrna), mrna_l95=lower(mrna), mrna_u95=upper(mrna),
              prot_mean=mean(prot), prot_l95=lower(prot), prot_u95=upper(prot)) %>%
    ungroup() -> panel_c_smry

  xlims = ylims_c = c(0.008, 2.3)
  plot(NA, NA, xlim=xlims, ylim=ylims_c, log='xy', axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, at=rep(1:9,2)*rep(10^(-2:-1),each=9), labels=NA, lwd=0, lwd.ticks=1, tck=-0.02)
  axis(side=1, at=c(0.01,0.1,1), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=c(0.01,0.1,1), labels=c('1','10','100'), lwd=0, line=-0.5, cex.axis=0.7)
  mtext(side=1, line=1.6, text=expression(residual~italic(PRNP)), cex=0.7)
  axis(side=2, at=ylims_c, labels=NA, lwd.ticks=0)
  axis(side=2, at=rep(1:9,2)*rep(10^(-2:-1),each=9), labels=NA, lwd=0, lwd.ticks=1, tck=-0.02)
  axis(side=2, at=c(0.01,0.1,1), labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=2, at=c(0.01,0.1,1), labels=c('1','10','100'), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=1.7, text='residual PrP', cex=0.7)
  abline(h=1, lty=2)
  abline(v=1, lty=3)
  suppressWarnings(arrows(x0=panel_c_smry$mrna_l95, x1=panel_c_smry$mrna_u95, y0=panel_c_smry$prot_mean, code=3, angle=90, length=0.02, lwd=1.5, col=panel_c_smry$color))
  suppressWarnings(arrows(x0=panel_c_smry$mrna_mean, y0=panel_c_smry$prot_l95, y1=panel_c_smry$prot_u95, code=3, angle=90, length=0.02, lwd=1.5, col=panel_c_smry$color))
  par(xpd=T)
  text(x=panel_c_smry$mrna_mean*0.5, y=panel_c_smry$prot_mean*1.6, srt=-45, labels=paste0(panel_c_smry$dose_disp,' µg'), col=panel_c_smry$color, cex=0.6)
  par(xpd=F)
  points(panel_c_data$mrna, panel_c_data$prot, pch=21, col=panel_c_data$color, bg='#FFFFFF')
  mtext(side=3, line=0.5, text='mRNA vs. PrP', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.25, line=0.2); panel = panel + 1

  write_supp_table(panel_c_smry %>% select(-color), 'Summary of transgenic mouse brain Prnp mRNA and PrP protein residual by ICV dose, 28 days post-dose (Figure 6A-C).')


  ### D. mouse ICV duration, mRNA, 300 µg ####
  par(mar=c(3,4,3,1))

  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='mRNA', dose_mg %in% c(0,0.3)) %>%
    inner_join(mouse_dose_meta %>% select(dose_mg,dose_disp,color), by='dose_mg') -> panel_d_data

  panel_d_data %>%
    group_by(dose_mg, dose_disp, color, days_post_dose) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() -> panel_d_smry

  xlims = c(20,92)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, at=c(28,84), labels=NA, lwd=0, lwd.ticks=1)
  axis(side=1, at=c(28,84), labels=c(28,84), lwd=0, line=-0.5, cex.axis=0.7)
  mtext(side=1, line=1.6, text='day post dose', cex=0.7)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text=expression(residual~italic(PRNP)), cex=0.7)
  abline(h=1, lty=2)
  for (this_dose in unique(panel_d_smry$dose_mg)) {
    subs = panel_d_smry %>% filter(dose_mg==this_dose) %>% arrange(days_post_dose)
    points(subs$days_post_dose, subs$mean, type='l', col=subs$color[1], lwd=1.5)
  }
  suppressWarnings(arrows(x0=panel_d_smry$days_post_dose, y0=panel_d_smry$l95, y1=panel_d_smry$u95, code=3, angle=90, length=0.04, col=panel_d_smry$color, lwd=1.5))
  points(panel_d_smry$days_post_dose, panel_d_smry$mean, pch=19, col=panel_d_smry$color, cex=0.9)
  set.seed(1)
  points(jitter(panel_d_data$days_post_dose, amount=2), panel_d_data$residual, pch=21, col=panel_d_data$color, bg='#FFFFFF')
  legend('topright', c('0 µg','300 µg'), col=mouse_dose_meta$color[c(1,5)], pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, lwd=1.5, bty='n', cex=0.6)
  mtext(side=3, line=0.5, text='mRNA', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.15, line=0.2); panel = panel + 1

  write_supp_table(panel_d_smry %>% select(-color), 'Summary of transgenic mouse brain Prnp mRNA residual over time, aCSF vs. 300 µg (Figure 6D).')


  ### E. mouse ICV duration, protein, 300 µg ####
  par(mar=c(3,4,3,1))

  alny_data %>%
    filter(species=='transgenic mouse', tissue=='brain', analyte=='protein', dose_mg %in% c(0,0.3)) %>%
    inner_join(mouse_dose_meta %>% select(dose_mg,dose_disp,color), by='dose_mg') -> panel_e_data

  panel_e_data %>%
    group_by(dose_mg, dose_disp, color, days_post_dose) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() -> panel_e_smry

  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=1, at=c(28,84), labels=NA, lwd=0, lwd.ticks=1)
  axis(side=1, at=c(28,84), labels=c(28,84), lwd=0, line=-0.5, cex.axis=0.7)
  mtext(side=1, line=1.6, text='day post dose', cex=0.7)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text='residual PrP', cex=0.7)
  abline(h=1, lty=2)
  for (this_dose in unique(panel_e_smry$dose_mg)) {
    subs = panel_e_smry %>% filter(dose_mg==this_dose) %>% arrange(days_post_dose)
    points(subs$days_post_dose, subs$mean, type='l', col=subs$color[1], lwd=1.5)
  }
  suppressWarnings(arrows(x0=panel_e_smry$days_post_dose, y0=panel_e_smry$l95, y1=panel_e_smry$u95, code=3, angle=90, length=0.04, col=panel_e_smry$color, lwd=1.5))
  points(panel_e_smry$days_post_dose, panel_e_smry$mean, pch=19, col=panel_e_smry$color, cex=0.9)
  set.seed(1)
  points(jitter(panel_e_data$days_post_dose, amount=2), panel_e_data$residual, pch=21, col=panel_e_data$color, bg='#FFFFFF')
  legend('topright', c('0 µg','300 µg'), col=mouse_dose_meta$color[c(1,5)], pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, lwd=1.5, bty='n', cex=0.6)
  mtext(side=3, line=0.5, text='PrP', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.15, line=0.2); panel = panel + 1

  write_supp_table(panel_e_smry %>% select(-color), 'Summary of transgenic mouse brain PrP protein residual over time, aCSF vs. 300 µg (Figure 6E).')


  ### F. NHP CSF, 84 days ####
  par(mar=c(4,4,3,1))

  alny_data %>%
    filter(species=='NHP', tissue=='CSF', !is.na(residual)) %>%
    inner_join(nhp_dose_meta, by='dose_mg') -> panel_f_data

  panel_f_data %>%
    group_by(dose_mg, dose_disp, color, days_post_dose) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() -> panel_f_smry

  ylims_f = c(0, 1.5)
  xlims_f = c(-5,90)
  plot(NA, NA, xlim=xlims_f, ylim=ylims_f, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims_f, labels=NA, lwd.ticks=0)
  axis(side=1, at=c(0,15,29,57,84), labels=NA, lwd=0, lwd.ticks=1)
  axis(side=1, at=c(0,15,29,57,84), labels=c(0,15,29,57,84), lwd=0, line=-0.5, cex.axis=0.65)
  mtext(side=1, line=1.6, text='day post dose', cex=0.7)
  axis(side=2, at=ylims_f, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:6/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:6/4, labels=percent(0:6/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text='residual PrP', cex=0.7)
  abline(h=1, lty=2)
  for (this_dose in unique(panel_f_smry$dose_mg)) {
    subs = panel_f_smry %>% filter(dose_mg==this_dose) %>% arrange(days_post_dose)
    points(subs$days_post_dose, subs$mean, type='l', col=subs$color[1], lwd=1.5)
  }
  suppressWarnings(arrows(x0=panel_f_smry$days_post_dose, y0=panel_f_smry$l95, y1=panel_f_smry$u95, code=3, angle=90, length=0.03, col=panel_f_smry$color, lwd=1.5))
  points(panel_f_smry$days_post_dose, panel_f_smry$mean, pch=19, col=panel_f_smry$color, cex=0.9)
  set.seed(1)
  points(jitter(panel_f_data$days_post_dose, amount=1.5), panel_f_data$residual, pch=21, col=panel_f_data$color, bg='#FFFFFF', cex=0.7)
  legend('topleft', nhp_dose_meta$dose_disp, col=nhp_dose_meta$color, pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, lwd=1.5, bty='n', cex=0.6)
  mtext(side=3, line=0.5, text='CSF PrP', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.09, line=0.2); panel = panel + 1

  write_supp_table(panel_f_smry %>% select(-color), 'Summary of NHP CSF PrP protein residual over time (Figure 6F).')


  ### G. NHP brain tissue mRNA, 85 days (renormalized to aCSF mean) ####
  par(mar=c(4,4,3,1))

  alny_data %>%
    filter(species=='NHP', analyte=='mRNA', tissue %in% region_meta$tissue) %>%
    group_by(tissue) %>%
    mutate(residual = residual / mean(residual[dose_mg==0])) %>%
    ungroup() %>%
    inner_join(nhp_dose_meta, by='dose_mg') %>%
    inner_join(region_meta, by='tissue') %>%
    mutate(x = region_x + offset) -> panel_g_data

  panel_g_data %>%
    group_by(tissue, region_x, dose_mg, dose_disp, color, x) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() -> panel_g_smry

  nhp_tissue_barwidth = 0.1
  xlims_g = c(0.5,4.5)
  plot(NA, NA, xlim=xlims_g, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims_g, labels=NA, lwd.ticks=0)
  mtext(side=1, at=region_meta$region_x, text=region_meta$tissue_disp, cex=0.55, line=0.0, padj=1)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=2.6, text=expression(residual~italic(PRNP)), cex=0.7)
  abline(h=1, lty=2)
  segments(x0=panel_g_smry$x-nhp_tissue_barwidth, x1=panel_g_smry$x+nhp_tissue_barwidth, y0=panel_g_smry$mean, col=panel_g_smry$color, lwd=2)
  suppressWarnings(arrows(x0=panel_g_smry$x, y0=panel_g_smry$l95, y1=panel_g_smry$u95, code=3, angle=90, length=0.03, col=panel_g_smry$color, lwd=1.5))
  set.seed(1)
  points(jitter(panel_g_data$x, amount=0.05), panel_g_data$residual, pch=21, col=panel_g_data$color, bg='#FFFFFF', cex=0.7)
  legend('topright', nhp_dose_meta$dose_disp, col=nhp_dose_meta$color, pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, bty='n', cex=0.55)
  mtext(side=3, line=0.5, text='brain mRNA', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.09, line=0.2); panel = panel + 1

  write_supp_table(panel_g_smry %>% select(-color,-x), 'Summary of NHP brain tissue Prnp mRNA residual by region, renormalized to aCSF mean (Figure 6G).')


  ### H. NHP brain tissue protein, 85 days ####
  par(mar=c(4,2.5,3,1))

  alny_data %>%
    filter(species=='NHP', analyte=='protein', tissue %in% region_meta$tissue) %>%
    inner_join(nhp_dose_meta, by='dose_mg') %>%
    inner_join(region_meta, by='tissue') %>%
    mutate(x = region_x + offset) -> panel_h_data

  panel_h_data %>%
    group_by(tissue, region_x, dose_mg, dose_disp, color, x) %>%
    summarize(.groups='keep', n=n(), mean=mean(residual,na.rm=T), l95=lower(residual), u95=upper(residual)) %>%
    ungroup() -> panel_h_smry

  plot(NA, NA, xlim=xlims_g, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims_g, labels=NA, lwd.ticks=0)
  mtext(side=1, at=region_meta$region_x, text=region_meta$tissue_disp, cex=0.55, line=0.0, padj=1)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  axis(side=2, at=0:5/4, labels=NA, lwd=0, lwd.ticks=1)
  axis(side=2, at=0:5/4, labels=percent(0:5/4), lwd=0, line=-0.5, las=2, cex.axis=0.7)
  mtext(side=2, line=1.75, text='residual PrP', cex=0.7)
  abline(h=1, lty=2)
  segments(x0=panel_h_smry$x-nhp_tissue_barwidth, x1=panel_h_smry$x+nhp_tissue_barwidth, y0=panel_h_smry$mean, col=panel_h_smry$color, lwd=2)
  suppressWarnings(arrows(x0=panel_h_smry$x, y0=panel_h_smry$l95, y1=panel_h_smry$u95, code=3, angle=90, length=0.03, col=panel_h_smry$color, lwd=1.5))
  set.seed(1)
  points(jitter(panel_h_data$x, amount=0.05), panel_h_data$residual, pch=21, col=panel_h_data$color, bg='#FFFFFF', cex=0.7)
  legend('topright', nhp_dose_meta$dose_disp, col=nhp_dose_meta$color, pch=21, pt.bg='#FFFFFF', pt.lwd=1.5, bty='n', cex=0.55)
  mtext(side=3, line=0.5, text='brain PrP', cex=0.7)
  mtext(LETTERS[panel], side=3, cex=1.3, adj=-0.09, line=0.2); panel = panel + 1

  write_supp_table(panel_h_smry %>% select(-color,-x), 'Summary of NHP brain tissue PrP protein residual by region (Figure 6H).')

  silence_is_golden = dev.off()
}


# SUPPLEMENT ####

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))


