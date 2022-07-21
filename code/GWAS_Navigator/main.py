# Todo redo the requirements!!

import pandas as pd

from data_provider import Model

from functools import partial

# bokeh stuff
from bokeh.layouts import column
from bokeh.plotting import figure, curdoc
from bokeh.models import LabelSet, ColorBar, Select, \
    CheckboxGroup, Span, Div, CustomJS, AutocompleteInput, \
    CrosshairTool, Button, RadioButtonGroup
from bokeh.transform import linear_cmap, factor_mark
from bokeh.palettes import Spectral6, Greys256
from bokeh.models import HoverTool, HBar, FactorRange, Title

'''Create to create the plots'''


def generate_manhatten_plot(data):
    # Create the plot
    manhatten_plot = figure(
        x_axis_label="Pos (on hg38)",
        y_axis_label="log₁₀(p)",
        x_axis_location="above",
        height=180,
        tools=['tap', 'reset', 'xwheel_zoom'],
        active_scroll="xwheel_zoom",
        y_axis_location="right",
    )

    # Add a line at p = 5*10^-8 (-> -log10(5*10^-8) = approx. 7.3)
    manhatten_plot.add_layout(Span(location=7.3,
                                   dimension='width', line_color='black',
                                   line_dash='dashed', line_width=1)
                              )

    # Scatter plot glyphes (-log p values)
    manhatten_plot.circle(x='Pos_hg38',
                          y='logP_value',
                          name='Scatter',
                          source=data.gwas_source,
                          size=6,
                          alpha=0.4,
                          color='gray',
                          )

    # Add overlay for consequences and LDs
    consequences = ['synonymous_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'upstream_gene_variant',
                    'downstream_gene_variant', 'missense_variant', 'stop_gained', 'start_lost', 'intron_variant',
                    'intergenic_variant', 'non_coding_transcript_exon_variant', 'non_coding_transcript_variant',
                    'splice_region_variant', 'splice_donor_variant', 'TF_binding_site_variant',
                    'regulatory_region_variant']
    marker_types = ['circle_dot', 'circle_x', 'circle_y', 'triangle', 'inverted_triangle', 'star', 'diamond_dot',
                    'diamond_cross', 'square_x', 'square_cross', 'square', 'square_dot', 'hex', 'hex_dot', 'star_dot',
                    'plus']

    # overlay glyph:
    ld_mapper = linear_cmap(field_name='r2', palette=Spectral6, low=0, high=1.0, nan_color='gray')
    overlay = manhatten_plot.scatter('Pos_hg38',
                                     'logP_value',
                                     name='Scatter',
                                     source=data.gwas_overlay_source,
                                     size=8,
                                     alpha=0.8,
                                     muted_alpha=0.4,
                                     color='black',
                                     line_alpha=1,
                                     line_width=0.3,
                                     fill_color=ld_mapper,
                                     marker=factor_mark('most_severe_consequence', markers=marker_types,
                                                        factors=consequences),
                                     )

    # Add hover tools:
    tooltip = """
    <div class="plot-tooltip">
        <div>
            <h4>@rsID</h4>
        </div>
        <div>
            <div>
                <span style="font-weight: bold;">Position: </span>@Pos_hg38
            </div>
            <div>
                <span style="font-weight: bold;">log₁₀(p): </span>@logP_value
            </div>
            <div>
                <span style="font-weight: bold;">MAF: </span>@maf
            </div>
            <div>
                <span style="font-weight: bold;">R2: </span>@r2
            </div>
            <div>
                <span style="font-weight: bold;">MSC: </span>@most_severe_consequence
            </div>
        </div>
    </div>
    """
    manhatten_plot.add_tools(HoverTool(renderers=[overlay], tooltips=tooltip))

    # Add color bar as legend for LDs:
    color_bar = ColorBar(color_mapper=ld_mapper['transform'],
                         label_standoff=12,
                         width=10,
                         title='LDs in selected Pop',
                         )
    manhatten_plot.add_layout(color_bar, 'right')

    manhatten_plot.scatter("Pos_hg38",
                           'logP_value',
                           source=data.gwas_proxy_source,
                           size=14, alpha=1,
                           muted_alpha=0.8,
                           color="purple",
                           line_alpha=1,
                           line_width=1.2,
                           fill_color=None,
                           marker="circle_dot",
                           )

    # Make everything pretty
    manhatten_plot.yaxis.major_tick_line_color = 'gray'
    manhatten_plot.yaxis.minor_tick_line_color = 'gray'
    manhatten_plot.xaxis.major_tick_line_color = 'gray'
    manhatten_plot.xaxis.minor_tick_line_color = 'gray'

    manhatten_plot.add_layout(Title(text="Coronary Artery Disease"), 'left')
    manhatten_plot.add_layout(Title(text="GWAS for"), 'left')

    return manhatten_plot


def generate_gene_plot(data, ref_plot):
    gene_plot = figure(x_range=ref_plot.x_range,
                       sizing_mode="stretch_width",
                       height=100,
                       tools=['tap', 'xwheel_zoom', 'hover'],
                       active_scroll="xwheel_zoom",
                       y_axis_location="right"
                       )

    renderer = gene_plot.hbar(y="level",
                              right="end",
                              left='start',
                              height=0.2,
                              source=data.gene_source,
                              fill_color='gray',
                              line_color='gray',
                              )

    selected_elem = HBar(
        fill_color='red',
        line_color='red',
        line_alpha=0.9,
        fill_alpha=0.9,
    )
    nonselected_elem = HBar(
        fill_color='gray',
        line_color='gray',
        line_alpha=0.7,
        fill_alpha=0.7,
    )

    renderer.selection_glyph = selected_elem
    renderer.nonselection_glyph = nonselected_elem

    tooltip = """
        <div class="plot-tooltip">
            <div>
                <h4>@name</h4>
            </div>
            <div>
                <div>
                    <span style="font-weight: bold;">Biotype: </span>@biotype
                </div>
                <div>
                    <span style="font-weight: bold;">Strand: </span>@strand
                </div>
            </div>
        </div>
        """

    hover = gene_plot.select(dict(type=HoverTool))
    hover.tooltips = tooltip

    # Add labels
    gene_labels = LabelSet(
        x='label_pos',
        y='level',
        text='name',
        x_offset=-5,
        y_offset=0,
        source=data.gene_source,
        text_font_size='7.2pt',
        text_align='center',
        text_color='color')
    gene_plot.add_layout(gene_labels)

    # Make everything pretty
    gene_plot.xaxis.visible = False
    gene_plot.xgrid.visible = True
    gene_plot.ygrid.visible = False
    gene_plot.yaxis.major_tick_line_color = None
    gene_plot.yaxis.axis_line_color = None
    gene_plot.toolbar.logo = None
    gene_plot.toolbar_location = None
    gene_plot.y_range.flipped = True
    gene_plot.yaxis.major_label_text_font_size = "7pt"

    gene_plot.add_layout(Title(text="Genes"), 'left')

    return gene_plot


def generate_reg_plot(data, ref_plot):
    reg_plot = figure(x_range=ref_plot.x_range,
                      sizing_mode="stretch_width",
                      height=10 * 7,
                      active_scroll="xwheel_zoom",
                      tools=['xwheel_zoom', 'hover'],
                      y_axis_location="right",
                      )

    reg_glyph = HBar(y="level",
                     right="end_GRCh38",
                     left='start_GRCh38',
                     height=0.7,
                     fill_color='color',
                     line_color=None,
                     fill_alpha=0.9,
                     )

    reg_plot.add_glyph(data.reg_build_source, reg_glyph)

    tooltip = """
        <div class="plot-tooltip">
            <div>
                <h4>@gene_id</h4>
            </div>
            <div>
                <div>
                    <span style="font-weight: bold;">Description: </span>@feature_type
                </div>
                <div>
                    <span style="font-weight: bold;">Strand: </span>@strand
                </div>
            </div>
        </div>
        """

    hover = reg_plot.select(dict(type=HoverTool))
    hover.tooltips = tooltip

    # Create labels and add them
    label_dict = {
        0: "Open chromatin region",
        1: "CTCF binding site",
        2: "Enhancer",
        3: "TSS",
        4: "TF binding site",
        5: "Promoter + Flanks",
        6: 'Constraints'
    }

    # ToDO: Add the TSS

    reg_plot.yaxis.ticker = [0, 1, 2, 3, 4, 5, 6]
    reg_plot.yaxis.ticker.num_minor_ticks = 0
    reg_plot.yaxis.major_label_overrides = label_dict

    # Make everything pretty
    reg_plot.yaxis.major_label_text_font_size = "7pt"
    reg_plot.yaxis.major_tick_line_color = None
    reg_plot.yaxis.axis_line_color = None
    reg_plot.xaxis.visible = False
    reg_plot.yaxis.minor_tick_line_color = None
    reg_plot.xgrid.visible = True
    reg_plot.toolbar.logo = None
    reg_plot.toolbar_location = None
    reg_plot.y_range.end = -0.5
    reg_plot.y_range.start = 6.5

    reg_plot.add_layout(Title(text="Elements"), 'left')
    reg_plot.add_layout(Title(text="reg."), 'left')

    return reg_plot


def generate_scATAC_seq_plot(data, ref_plot):
    seq_plot = figure(
        x_range=ref_plot.x_range,
        y_range=FactorRange(factors=[]),
        sizing_mode="stretch_width",
        height=100,
        active_scroll="xwheel_zoom",
        tools=['xwheel_zoom', 'hover'],
        y_axis_location="right"
    )

    # Color mapper (gray values) for the score
    score_mapper = linear_cmap(field_name='score', palette=Greys256, low=1000, high=0)

    bar_glyph = HBar(y='biosample',
                     right="end_hg38",
                     left='start_hg38',
                     height=0.4,
                     fill_color=score_mapper,
                     line_color='gray',
                     line_alpha=0.5,
                     fill_alpha=0.9,
                     )

    summit_glyph = HBar(y='biosample',
                        right='center',
                        left='center',
                        line_color='red',
                        line_width=0.2,
                        line_alpha=0.8,
                        height=0.4
                        )

    seq_plot.add_glyph(data.Miller_etal_source, bar_glyph)
    seq_plot.add_glyph(data.Miller_etal_source, summit_glyph)

    tooltip = """
            <div class="plot-tooltip">
                <div>
                    <h4>@peak_name</h4>
                </div>
                <div>
                    <div>
                        <span style="font-weight: bold;">Target: </span>@target
                    </div>
                    <div>
                        <span style="font-weight: bold;">Score: </span>@score
                    </div>
                </div>
            </div>
            """

    hover = seq_plot.select(dict(type=HoverTool))
    hover.tooltips = tooltip

    # Make everything pretty
    seq_plot.yaxis.major_label_text_font_size = "6pt"
    seq_plot.yaxis.major_tick_line_color = None
    seq_plot.yaxis.axis_line_color = None
    seq_plot.xaxis.visible = False
    seq_plot.yaxis.minor_tick_line_color = None
    seq_plot.xgrid.visible = True
    seq_plot.toolbar.logo = None
    seq_plot.toolbar_location = None

    seq_plot.add_layout(Title(text="seq"), 'left')
    seq_plot.add_layout(Title(text="scATAC"), 'left')

    return seq_plot


def generate_catlas_plot(data, ref_plot):
    seq_plot = figure(
        x_range=ref_plot.x_range,
        y_range=FactorRange(factors=[]),
        sizing_mode="stretch_width",
        height=100,
        active_scroll="xwheel_zoom",
        tools=['xwheel_zoom', 'hover'],
        y_axis_location="right"
    )

    # Color mapper (gray values) for the score
    score_mapper = linear_cmap(field_name='score', palette=Greys256, low=1000, high=0)

    bar_glyph = HBar(y='biosample',
                     right="end_hg38",
                     left='start_hg38',
                     height=0.4,
                     fill_color=score_mapper,
                     line_color='gray',
                     line_alpha=0.5,
                     fill_alpha=0.9,
                     )

    summit_glyph = HBar(y='biosample',
                        right='center',
                        left='center',
                        line_color='red',
                        line_width=0.2,
                        line_alpha=0.8,
                        height=0.4
                        )

    seq_plot.add_glyph(data.CATlas_source, bar_glyph)
    seq_plot.add_glyph(data.CATlas_source, summit_glyph)

    TOOLTIP = """
            <div class="plot-tooltip">
                <div>
                    <h4>@biosample</h4>
                </div>
                <div>
                    <div>
                        <span style="font-weight: bold;">Score: </span>@score
                    </div>
                </div>
            </div>
            """

    hover = seq_plot.select(dict(type=HoverTool))
    hover.tooltips = TOOLTIP

    # Make everything pretty
    seq_plot.yaxis.major_label_text_font_size = "6pt"
    seq_plot.yaxis.major_tick_line_color = None
    seq_plot.yaxis.axis_line_color = None
    seq_plot.xaxis.visible = False
    seq_plot.yaxis.minor_tick_line_color = None
    seq_plot.xgrid.visible = True
    seq_plot.toolbar.logo = None
    seq_plot.toolbar_location = None

    seq_plot.add_layout(Title(text="CATlas"), 'left')

    return seq_plot


def generate_tads_plot(data, ref_plot):
    tads_plot = figure(x_range=ref_plot.x_range,
                       y_range=FactorRange(factors=[]),
                       sizing_mode="stretch_width",
                       height=100,  # inital heights needs to be fixed !!!
                       active_scroll="xwheel_zoom",
                       tools=['xwheel_zoom', 'hover'],
                       y_axis_location="right"
                       )

    bar_glyph = HBar(y='biosample',
                     right='start_hg38',
                     left='end_hg38',
                     height=0.4,
                     fill_color='gray',
                     line_alpha=0.5,
                     fill_alpha=0.9,
                     )

    tads_plot.add_glyph(data.TADs_source, bar_glyph)

    TOOLTIP = """
            <div class="plot-tooltip">
                <div>
                    <h6>@biosample</h6>
                </div>
            </div>
            """

    hover = tads_plot.select(dict(type=HoverTool))
    hover.tooltips = TOOLTIP

    # Make everything pretty
    tads_plot.yaxis.major_label_text_font_size = "6pt"
    tads_plot.yaxis.major_tick_line_color = None
    tads_plot.yaxis.axis_line_color = None
    tads_plot.xaxis.visible = False
    tads_plot.yaxis.minor_tick_line_color = None
    tads_plot.xgrid.visible = True
    tads_plot.toolbar.logo = None
    tads_plot.toolbar_location = None

    tads_plot.add_layout(Title(text="TADs"), 'left')

    return tads_plot


def generate_trait_plot(data, ref_plot):
    trait_plot = figure(x_range=ref_plot.x_range,
                        sizing_mode="stretch_width",
                        y_range=FactorRange(factors=[]),
                        height=100,
                        active_scroll="xwheel_zoom",
                        tools=['xwheel_zoom', 'hover'],
                        y_axis_location="right"
                        )

    trait_plot.scatter(x='pos_Hg38',
                       y='trait',
                       name='Scatter',
                       source=data.traits_source,
                       size='size',
                       alpha=1,
                       color='color',
                       marker='marker'
                       )

    tooltip = f"""
            <div class="plot-tooltip">
                <div>
                    <h4>@trait</h4>
                </div>
                <div>
                    <div>
                        <span style="font-weight: bold;">rsID: </span>rs@rsid
                    </div>
                    <div>
                        <span style="font-weight: bold;">pos on Hg38: </span>@chr_Hg38:@pos_Hg38
                    </div>
                    <div>
                        <span style="font-weight: bold;">p-value Exponent: </span>@p_value
                    </div>
                    <div>
                        <span style="font-weight: bold;">Beta: </span>@beta
                    </div>
                    <div>
                        <span style="font-weight: bold;">Reported Gene(s): </span>@reported_genes
                    </div>

                    
                </div>
            </div>
            """

    hover = trait_plot.select(dict(type=HoverTool))
    hover.tooltips = tooltip

    # Make everything pretty
    trait_plot.yaxis.major_label_text_font_size = "6pt"
    trait_plot.yaxis.major_tick_line_color = None
    trait_plot.yaxis.axis_line_color = None
    trait_plot.xaxis.visible = False
    trait_plot.yaxis.minor_tick_line_color = None
    trait_plot.xgrid.visible = True
    trait_plot.toolbar.logo = None
    trait_plot.toolbar_location = None

    trait_plot.add_layout(Title(text="Traits"), 'left')

    return trait_plot


def generate_abc_plot(data, ref_plot):
    abc_plot = figure(
        x_range=ref_plot.x_range,
        y_range=FactorRange(factors=[]),
        sizing_mode="stretch_width",
        height=100,
        active_scroll="xwheel_zoom",
        tools=[
            'xwheel_zoom',
            'hover'
        ],
        y_axis_location="right"
    )

    # Color mapper (gray values) for the score
    element = abc_plot.hbar(y='CellType',
                            right="end_hg38",
                            left='start_hg38',
                            height=0.5,
                            source=data.ABC_source,
                            fill_color='gray',
                            line_color='gray',
                            )

    selected_elem = HBar(
        fill_color='red',
        line_color='red',
        line_alpha=0.9,
        fill_alpha=0.9,
    )
    nonselected_elem = HBar(
        fill_color='gray',
        line_color='gray',
        line_alpha=0.7,
        fill_alpha=0.7,
    )

    element.selection_glyph = selected_elem
    element.nonselection_glyph = nonselected_elem

    tooltip = """
        <div class="plot-tooltip">
            <div>
                <h4>@CellType</h4>
            </div>
            <div>
                <div>
                    <span style="font-weight: bold;">Target Gene: </span>@TargetGene
                </div>
                <div>
                    <span style="font-weight: bold;">ABC Score: </span>@Score
                </div>
            </div>
        </div>
        """

    hover = abc_plot.select(dict(type=HoverTool))
    hover.tooltips = tooltip

    abc_plot.yaxis.major_label_text_font_size = "6pt"
    abc_plot.yaxis.major_tick_line_color = None
    abc_plot.yaxis.axis_line_color = None
    abc_plot.xaxis.visible = False
    abc_plot.yaxis.minor_tick_line_color = None
    abc_plot.xgrid.visible = True
    abc_plot.toolbar.logo = None
    abc_plot.toolbar_location = None

    abc_plot.add_layout(Title(text="ABC"), 'left')

    return abc_plot


'''Create the update functions for stuff happening in the plot'''


# These functions are important to trigger an update on the feedback div to display a text
# that is saved in the settings manager
def update_feedback_div(model):
    feedback_div.text = model.feedback_text


def update_gene_div(model):
    gene_div.text = model.gene_text


def snp_input_event(attr, old, new):
    if new != old:
        model.update_feedback_div_text(f'⌛ Fetching info for {new}')
        update_feedback_div(model)

        # schedule the callback which updates the parameters in the settings Manager
        # https://discourse.bokeh.org/t/show-loading-sign-during-calculations/4410
        model.update_snp(new)
        curdoc().add_next_tick_callback(partial(update_feedback_div, model))

        # Update the text of the button to query.
        # This needs to be done in the next iteration of the eventloop when the updates have happend!
        def update_the_snpbutton():
            if model.rsid == '':
                snp_button.disabled = True
                snp_button.label = 'Please enter a valid rsID.'
            else:
                snp_button.disabled = False
                snp_button.label = f'Generate plot for {model.rsid}'

        curdoc().add_next_tick_callback(update_the_snpbutton)


def gene_input_event(attr, old, new):
    if new != old:
        model.update_gene_div_text(f'⌛ Fetching info for {new}')
        update_gene_div(model)

        # schedule the callback which updates the parameters in the settings Manager
        # https://discourse.bokeh.org/t/show-loading-sign-during-calculations/4410
        curdoc().add_next_tick_callback(partial(model.update_gene, new))
        curdoc().add_next_tick_callback(partial(update_gene_div, model))

        # Update the buttons to query from gene search
        # This needs to be done in the next iteration of the eventloop when the updates have happend!
        def update_the_genebutton():
            if len(model.close_snps) == 0:
                gene_button.disabled = True
                gene_button.labels = ['Please enter a viable gene name.']
            else:
                gene_button.disabled = False
                gene_button.labels = model.close_snps.tolist()
                gene_button.active = None

        curdoc().add_next_tick_callback(update_the_genebutton)


def update_plot_from_snpinput(data_provider, plotlist):
    model.update_feedback_div_text(f'⌛ Loading data for {model.rsid}, this might take a while.')
    update_feedback_div(model)


    # General function to UPDATE THE PLOT (reused for input of a gene)
    curdoc().add_next_tick_callback(partial(update_plot, data_provider, plotlist))

    # Update the div
    curdoc().add_next_tick_callback(partial(update_feedback_div, model))

    # Reset the gene query button
    gene_button.active = None
    snp_button.disabled = True


def update_plot_from_geneinput(active, data_provider, plotlist):
    if active is not None:
        model.update_snp(model.close_snps[active])

        if model.rsid == '':
            model.update_gene_div_text(f'WTF, it seems like you selected a variant that is not part of the GWAS summary statistics. How did we even get here (Probably because some versions are not synchronized and this exact variant got renamed...)?! I am very sorry.')
            curdoc().add_next_tick_callback(partial(update_feedback_div, model))
            curdoc().add_next_tick_callback(partial(update_gene_div, model))
            return

        # Update the text of the gene_div & then update the div in the plot
        model.update_gene_div_text(f'⌛ Loading data for {model.rsid}, this might take a while.')
        update_gene_div(model)

        # General function to UPDATE THE PLOT (reused for input of a SNP)
        curdoc().add_next_tick_callback(partial(update_plot, data_provider, plotlist))
        curdoc().add_next_tick_callback(partial(update_feedback_div, model))
        curdoc().add_next_tick_callback(partial(update_gene_div, model))

        # And update the SNP input
        snp_input.value = model.close_snps[active]


def update_plot(data_provider, plotlist):
    # Update the data
    data_provider.update_all_data()

    # traits
    traits = data_provider.traits_source.data['trait'].unique().tolist()
    if not traits:
        plotlist['trait_plot'].visible = False
    else:
        plotlist['trait_plot'].visible = True
        plotlist['trait_plot'].y_range.factors = traits
        plotlist['trait_plot'].height = len(traits) * 10

    # genes
    label_dict = pd.DataFrame(data_provider.gene_source.data).groupby(
        'biotype',
        as_index=False
    )['level'].min().set_index('level')['biotype'].to_dict()
    max_level = data_provider.gene_source.data['level'].max()

    plotlist['gene_plot'].yaxis.ticker = list(label_dict.keys())
    plotlist['gene_plot'].yaxis.major_label_overrides = label_dict
    plotlist['gene_plot'].y_range.end = -1
    plotlist['gene_plot'].y_range.start = max_level + 1

    # reg elm should not require an update

    # TADs
    tads = data_provider.TADs_source.data['biosample'].unique().tolist()
    if not tads:
        plotlist['tad_plot'].visible = False
    else:
        plotlist['tad_plot'].visible = True
        plotlist['tad_plot'].y_range.factors = tads
        plotlist['tad_plot'].height = len(tads) * 8

    # catlas
    biosamples = data_provider.CATlas_source.data['biosample'].unique().tolist()
    if not biosamples:
        plotlist['catlas_plot'].visible = False
    else:
        plotlist['catlas_plot'].visible = True
        plotlist['catlas_plot'].y_range.factors = biosamples
        plotlist['catlas_plot'].height = len(biosamples) * 10

    # scATAC
    biosamples = data_provider.Miller_etal_source.data['biosample'].unique().tolist()
    if not biosamples:
        plotlist['scatac_plot'].visible = False
    else:
        plotlist['scatac_plot'].visible = True
        plotlist['scatac_plot'].y_range.factors = biosamples
        plotlist['scatac_plot'].height = len(biosamples) * 10

    # ABC
    biosamples = data_provider.ABC_source.data['CellType'].unique().tolist()
    if not biosamples:
        plotlist['abc_plot'].visible = False
    else:
        plotlist['abc_plot'].visible = True
        plotlist['abc_plot'].y_range.factors = biosamples
        plotlist['abc_plot'].height = len(biosamples) * 10


'''Create the update functions for the settings tab'''


def pop_select_event(attr, old, new):
    if new != old:
        model.update_population(new)


def update_catlas_plot(active, data_provider, plotlist):
    all_samples_dict = {idx: sample for idx, sample in enumerate(data_provider.all_catlas_biosample)}

    active_samples = [all_samples_dict[idx] for idx in active]

    if not active_samples:
        plotlist['catlas_plot'].visible = False

    else:
        plotlist['catlas_plot'].visible = True
        plotlist['catlas_plot'].y_range.factors = active_samples
        plotlist['catlas_plot'].height = len(active_samples) * 10

    return


def update_atac_seq_plot(active, data_provider, plotlist):
    all_samples_dict = {idx: sample for idx, sample in enumerate(data_provider.all_miller_biosample)}

    active_samples = [all_samples_dict[idx] for idx in active]

    if not active_samples:
        plotlist['scatac_plot'].visible = False

    else:
        plotlist['scatac_plot'].visible = True
        plotlist['scatac_plot'].y_range.factors = active_samples
        plotlist['scatac_plot'].height = len(active_samples) * 10

    return


def update_abc_plot(active, data_provider, plotlist):
    all_samples_dict = {idx: sample for idx, sample in enumerate(data_provider.all_abc_celltypes)}

    active_samples = [all_samples_dict[idx] for idx in active]
    if not active_samples:
        plotlist['abc_plot'].visible = False

    else:
        plotlist['abc_plot'].visible = True
        plotlist['abc_plot'].y_range.factors = active_samples
        plotlist['abc_plot'].height = len(active_samples) * 10

    return


'''Create the update functions for the tap tool'''
# Todo: implement


'''Create the data Provider'''

model = Model()

'''Create the plots'''
manhatten_plot = generate_manhatten_plot(model)
gene_plot = generate_gene_plot(model, manhatten_plot)
reg_plot = generate_reg_plot(model, manhatten_plot)
scatac_plot = generate_scATAC_seq_plot(model, manhatten_plot)
catlas_plot = generate_catlas_plot(model, manhatten_plot)
tad_plot = generate_tads_plot(model, manhatten_plot)
trait_plot = generate_trait_plot(model, manhatten_plot)
abc_plot = generate_abc_plot(model, manhatten_plot)

plotlist = {
    'manhatten_plot': manhatten_plot,
    'trait_plot': trait_plot,
    'gene_plot': gene_plot,
    'reg_plot': reg_plot,
    'tad_plot': tad_plot,
    'scatac_plot': scatac_plot,
    'catlas_plot': catlas_plot,
    'abc_plot': abc_plot
}

# Add the crosshair to all plots:
crosshair = CrosshairTool(dimensions="height")
for plot in plotlist.values():
    plot.add_tools(crosshair)

# Add a reset to the x-range after change of the displayed SNPs.
# I think me scheduling the updates to the event loop to include update to the feedback divs
# messes with the normal updates.
reset_callback = CustomJS(args=dict(p=manhatten_plot), code="""p.reset.emit()""")
model.gwas_source.js_on_change('data', reset_callback)

'''Create all the input widgets:'''
feedback_div = Div(
    text=model.feedback_text,
    name="feedback_div",
    sizing_mode="stretch_width"
)

gene_div = Div(
    text=model.gene_text,
    name='gene_div',
    sizing_mode="stretch_width"
)

snp_input = AutocompleteInput(
    title='Select a SNP',
    name='snp_input',
    case_sensitive=False,
    min_characters=1,
    restrict=False
)
snp_input.on_change('value', snp_input_event)

gene_input = AutocompleteInput(
    title='Enter a gene',
    name='gene_input',
    case_sensitive=False,
    min_characters=1,
    restrict=False
)
gene_input.on_change('value', gene_input_event)

snp_button = Button(
    label="Generate plot",
    name='snp_button',
    button_type="default",
    disabled=True
)
snp_button.on_click(partial(update_plot_from_snpinput, model, plotlist))

gene_button = RadioButtonGroup(
    labels=['please search for a gene first.'],
    name='update_gene_rbutton',
    button_type="default",
    disabled=True
)
gene_button.on_click(partial(update_plot_from_geneinput, data_provider=model, plotlist=plotlist))

'''Create all the controls in the settings tab'''
population_select = Select(
    title="Population for LD structure:",
    value=model.population,
    options=['ALL', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
    name='pop_select',
    # disabled=True
)
population_select.on_change('value', pop_select_event)

atac_select = CheckboxGroup(
    labels=model.all_miller_biosample,
    name='atac_seq_select',
    # disabled=True,
    active=list(range(0, len(model.all_miller_biosample)))
)
atac_select.on_click(partial(update_atac_seq_plot, data_provider=model, plotlist=plotlist))

catlas_select = CheckboxGroup(
    labels=model.all_catlas_biosample,
    name='catlas_select',
    # disabled=True,
    active=list(range(0, len(model.all_catlas_biosample)))
)
catlas_select.on_click(partial(update_catlas_plot, data_provider=model, plotlist=plotlist))

abc_select = CheckboxGroup(
    labels=model.all_abc_celltypes,
    name='abc_select',
    # disabled=True,
    active=list(range(0, len(model.all_abc_celltypes)))
)
abc_select.on_click(partial(update_abc_plot, data_provider=model, plotlist=plotlist))

'''Add everything to the document'''
curdoc().add_root(column(children=list(plotlist.values()), name='plots', sizing_mode='scale_width'))

curdoc().add_root(feedback_div)
curdoc().add_root(gene_div)
curdoc().add_root(snp_input)
curdoc().add_root(gene_input)
curdoc().add_root(snp_button)
curdoc().add_root(gene_button)

curdoc().add_root(population_select)
curdoc().add_root(abc_select)
curdoc().add_root(atac_select)
curdoc().add_root(catlas_select)

# Initialize the plot with a SNP
model.update_snp('rs1230666')
update_plot(model, plotlist)