
draw_overall_menu = function(clusters) {
    complete_with_header(details_data);
    var upper = d3.select("#upper").append('div')
        .attr("class", "content");

    // title
    var content = upper.append('div')
        .attr('class', 'content')
        .html('<h2 >antiSMASH run ID</h2>'  +
              '<button id="btn-direction" type="button" class="btn btn-danger btn-xs" data-toggle="tooltip" title="Swap ORFs position for minus strand">Swap ORFs</button>' +
              '<p>'+ run_id + '</p>');

    // clusters button
    content.selectAll('.cluster-buttn')
        .data(clusters.sort(d3.ascending))
        .enter()
        .append('button')
        .attr('class', 'cluster-buttn')
        .text(function(d) {
            return d;
        })
        .on("click", function(d, index) {
            window.active_cluster = d;
            update_view_on_cluster_selection(d, this);
        });
}


cluster_draw_svgs = function(cluster, height, max_width) {
    margin = 20;
    domain_height = 40
    label_font_size = 14;
    svg_min_width = 65;

    d3.select("#cluster-details").selectAll('svg').remove();
    var container = d3.select("#cluster-details")
        .attr('class', 'content');

    // Determine xscale
    max_orf_len = 0;
    for (i=0; i < cluster['orfs'].length; i++) {
        max_orf_len = Math.max(max_orf_len, cluster['orfs'][i].sequence.length);
    }

    var x = d3.scale.linear()
        .domain([1, max_orf_len])
        .range([0, max_width]);

    // Draw one SVG by ORF
    var svgs = container.selectAll('svg')
        .data(cluster['orfs'])
        .enter()
        .append('svg')
        .attr("height", height)
        .attr("width", function (orf) { return Math.max(svg_min_width, x(orf.sequence.length + margin)) })
        .attr("id", function(d) { return d['id'] })
        .attr("class", "orf")
        .on("mousedown", function(d, index) { window.active_orf = d['id'] });

    // Draw the upper left label & the horizontal line
    svgs.append('text')
        .text(function (orf) { return orf['id'] })
        .attr("x", 5)
        .attr("y", label_font_size/2 + 5)
        .attr("class", "orf-label");

    svgs.append('line')
        .attr("stroke", function(orf) { return ((get_orf_strand(geneclusters, cluster['id'], orf['id']) < 0) ? 'red': 'black') })
        .attr("x1", 0)
        .attr("y1", margin + (domain_height/2))
        .attr("x2", function (orf) { return x(orf.sequence.length) })
        .attr("y2", margin + (domain_height/2))
        .attr("class", "orf-line");

    // Draw domains (rect, then label)
    var domains_group = svgs.selectAll('g')
        .data(function(orf) { return orf['domains'] })
        .enter()
        .append('g')
        .attr('class', 'domain')
        .on("click", function(d,i) { update_view_on_domain_selection(d,i, this) } );

    var domains_rect = domains_group
        .append("rect")
        .attr("x", function(d) { return x(d.start) })
        .attr("y", margin)
        .attr("rx", 14)
        .attr("ry", 14)
        .attr("width", function(d) { return x(d.end) - x(d.start) })
        .attr("height", domain_height)
        .attr("id", function(d, idx) { return "domain-rect" + idx })
        .attr("class", "domain-rect")
        .attr("fill", function(d) { return domain_fill_color(d['type']);  })
        .attr("stroke", function(d) { return domain_stroke_color(d['type']);  })
        .attr("stroke-width", 1);

    var domains_label_type = domains_group
        .append("text")
        .text(function(d) { return domain_label_type(d['type']); })
        .attr("x", function(d) { return (x(d.start) + (x(d.end) - x(d.start)) / 2) })
        .attr("y", margin + domain_height/2)
        .attr("text-anchor", "middle") //centering text horizontally
        .attr("dominant-baseline", "central") //centering vertically
    //Update font-size to stay into the rect boundaries
        .style("font-size", function(domain) { return Math.min(label_font_size, (x(domain.sequence.length / this.getComputedTextLength()) * label_font_size)) + "px"; })
        .attr("id", function(d, idx) { return "domain-label-type-" + idx})
        .attr("class", "domain-label-type")
    ;//.attr("font-weight", "bold");

    var domains_label_pred = domains_group
        .append("text")
        .text(function(d) { return domain_prediction(d) })
        .attr("x", function(d) { return x(d.start) + (x(d.end) - x(d.start))/2})
        .attr("text-anchor", "middle") //centering text horizontally
        .attr("dominant-baseline", "central") //centering vertically
        .attr("y", margin + domain_height + margin)
        .attr("class", "domain-label-pred");
}


update_view_on_cluster_selection = function(cluster, button) {
    // remove all content
    d3.select('#lower').selectAll('.content').remove();
    // deselect cluster button
    d3.selectAll('.cluster-buttn').attr('class', 'cluster-buttn');
    d3.select(button).attr('class', 'cluster-buttn active');
    // draw
    cluster_draw_summary(geneclusters[cluster]);
    cluster_draw_svgs(details_data[cluster], 100, 700);
}


update_view_on_domain_selection = function(domain, idx, element) {

    // deselect domain
    d3.selectAll('g.domain').attr('class', 'domain');
    d3.select(element).attr('class', 'domain active');

    domain_draw_summary(domain, idx)
}


retrieve_sequence_name = function(cluster_name, geneclusters) {
    /*
     *  Retrieve the sequence name of the sequence used to predict
     *  the cluster
     *  As this field is optional, return undef if we can't retrieve it
     */
    var seq_pattern = /id=(.+?)&/;
    var seq_match = seq_pattern.exec(geneclusters[cluster_name]['orfs'][0]['description']);
    var seq_name = seq_match[1];
    return seq_name;
}

cluster_draw_summary = function(cluster) {

    d3.select("#cluster-summary").selectAll('div').remove();

    var container = d3.select("#cluster-summary");
    container.attr('class', 'content');

    var panel = container.append('div')
        .attr('class', 'panel panel-primary');

    panel.append('div')
        .attr('class', 'panel-heading')
        .append('h3')
        .attr('class', 'panel-title')
        .text('Cluster-' + cluster['idx']);

    var seq_name = retrieve_sequence_name('cluster-' + cluster['idx'], geneclusters);
    panel.append('div')
        .attr('class', 'panel-body')
        .html('Type: ' + cluster['type'] + '<br>' +
              'Probability: ' + cluster['probability'] + '<br>' +
              'Sequence: ' + seq_name + '<br>' +
              'Location: ' + cluster['start'] + '-' + cluster['end']
             )
}


domain_prediction = function(domain) {
    var label = domain_label_type(domain['type']);
    var predictions =  domain['predictions'];


    switch (label) {
    case "A":
        return predictions[predictions.length - 1][1]
    case "C":
        return domain['type'].split('_')[1]
    default:
        return ''
    }
}


domain_draw_summary = function(domain, idx) {
    window.active_domain = 'domain_' + idx;

    var lower = d3.select('#lower');
    lower.selectAll('div').remove();


    var container = lower.append('div')
        .attr('class', 'content');

    // draw napdos results
    if ('napdos_details' in domain){
        container.append('h2')
            .text('NaPDoS table');

        var table = container.append('table')
            .attr('class', 'table table-bordered');
        var thead = table.append('thead');
        var tbody = table.append('tbody');

        var head_tr = thead.selectAll()
            .data(domain['napdos_details'])
            .enter()
            .append('tr')
            .filter(function (d, i) { return i === 0;});

        var body_tr = tbody.selectAll()
            .data(domain['napdos_details'])
            .enter()
            .append('tr')
            .filter(function (d, i) { return i != 0;});

        // Now create the table cells
        head_tr.selectAll("td")
            .data(function(d) {return d; })
            .enter()
            .append("th")
            .attr('class', 'word-wrap')
            .text(function(d) {return d;});
        body_tr.selectAll("td")
            .data(function(d) {return d; })
            .enter()
            .append("td")
            .text(function(d) {return d;});

        // Add href to NaPDoS svg tree
        container.append('a')
            .attr('href', 'napdos/svg/' + active_cluster)
            .attr('target', '_blank')
            .text('NaPDoS Tree');

    }

    // draw nrps results
    if ('nrpspredictor_details' in domain){
	var half_div = container.append('div')
            .attr('class', 'half');

	half_div.append('h2')
            .text('NRPSPredictor table');
        console.log(domain['nrpspredictor_details']);
	for (var i = 0, len = domain['nrpspredictor_details'].length; i < len; i++) {
	    if (i === 0) continue;
	    var row = 	domain['nrpspredictor_details'][i];
	    var table = half_div.append('table')
		.attr('class', 'table table-bordered');
            var tbody = table.append('tbody');
	    //.attr('class', "text-center")

	    var outside_applicable_domain = row[10];
	    tbody.append('tr').append('th').attr('colspan', 2).html(function(value) {
		var html = row[0];
		if (outside_applicable_domain == 1) {
		    html += '<span class="glyphicon glyphicon-remove"></span>';
		}
		else {
		    html += '<span class="glyphicon glyphicon-ok"></span>';
		}
		return html;
	    });
	    var tr  = tbody.append('tr');
	    tr.append('th').text("Signature")
	    tr.append('td').text( row[1] + ' / ' + row[2]);
	    tr  = tbody.append('tr').attr('bgcolor', "#E8E8E8");
	    tr.append('th').append('u').text("NRPSpredictor1");
	    tr.append('th').append('u').text("Predictions");
	    tr  = tbody.append('tr');
	    tr.append('th').text('Large Clusters');
	    tr.append('td').text(row[8]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('Small Clusters');
	    tr.append('td').text(row[9]);
	    tr  = tbody.append('tr').attr('bgcolor', "#E8E8E8");
	    tr.append('th').append('u').text("NRPSpredictor2");
	    tr.append('th').append('u').text("Predictions");
	    tr  = tbody.append('tr');
	    tr.append('th').text('Three Clusters');
	    tr.append('td').text(row[3]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('Large Clusters');
	    tr.append('td').text(row[4]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('Small Clusters');
	    tr.append('td').text(row[5]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('Single AA');
	    tr.append('td').text(row[6]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('Nearest Neighbor');
	    tr.append('td').text(row[7])
	    tr  = tbody.append('tr');
	    tr.append('th').text('Location');
	    tr.append('td').text(row[11]);
	    tr  = tbody.append('tr');
	    tr.append('th').text('PFAM score');
	    tr.append('td').text(row[12]);
	    tr  = tbody.append('tr');
	}

    }

    // add prediction table for A domains
    if ('predictions' in domain && domain['predictions'].length) {
        var half_div = container.append('div')
            .attr('class', 'half');

        half_div.append('h2')
            .text('Substrate Predictions:');

        var pred_tbody = half_div.append('div')
            .attr('class', 'half')
            .append('table')
            .attr('class', 'table table-bordered')
            .append('tbody');

        pred_tbody.selectAll('tr')
            .data(domain['predictions'])
            .enter()
            .append('tr')
            .html(function(d) { return '<th>' + d[0] + '</th>' + '<td>' + d[1] + '</td>' } );
    }

    // add fasta sequence
    var half_div = container.append('div')
        .attr('class', 'half');

    half_div.append('h2')
        .text('Domain\'s sequence');

    // copy to clipboard button
    var div = half_div.append('div').attr('class', "input-group");

    var fasta = domain['header'] + '\n'
        + chunk_string(domain['sequence'], 80);

    var span = div.append('span')
        .attr('id', "clipboard")
	.attr('class', "input-group-addon");

    span.append('i')
	.attr('class', "glyphicon glyphicon-copy glyphicon-align-left");

    span.on("click", function() {

        $("#sequence").select();
        document.execCommand('copy');
	//$("#sequence").blur();
        });

    div.append('textarea')
        .attr('id', 'sequence')
        .attr('readonly', '')
        .attr('rows', 4)
	.attr('class', "message form-control")
        .text(fasta);

    container.append('footer').attr('class', "footer");
}


domain_stroke_color = function(type) {
    switch (type) {
    case "AMP-binding":
    case "AOX":
        return "rgb(87,22,128)"
    case "PCP":
    case "ACP":
        return "rgb(11,78,199)"
    case "Cglyc":
    case "CXglyc":
    case "Condensation_DCL":
    case "Condensation_LCL":
    case "Condensation_Starter":
    case "Condensation_Dual":
    case "Heterocyclization":
        return "rgb(59,59,140)"
    case "Epimerization":
        return "rgb(59,59,140)"
    case "NRPS-COM_Nterm":
    case "NRPS-COM_Cterm":
    case "PKS_Docking_Nterm":
    case "PKS_Docking_Cterm":
    case "Trans-AT_docking":
        return "rgb(71,71,159)"
    case "Thioesterase":
    case "TD":
        return "rgb(119,3,116)"
    case "PKS_KS":
        return "rgb(9,179,9)"
    case "PKS_AT":
        return "rgb(221,6,6)"
    case "PKS_KR":
        return "rgb(10,160,76)"
    case "PKS_DH":
    case "PKS_DH2":
    case "PKS_DHt":
        return "rgb(186,103,15)"
    case "PKS_ER":
        return "rgb(12,161,137)"
    case "Aminotran_1_2":
    case "Aminotran_3":
    case "Aminotran_4":
    case "Aminotran_5":
    case "Polyketide_cyc2":
    default:
        return "rgb(147,147,147)"
    }
}


domain_fill_color = function(type) {
    switch (type) {
    case "AMP-binding":
    case "AOX":
        return "rgb(188,127,245)"
    case "PCP":
    case "ACP":
        return "rgb(129,190,247)"
    case "Cglyc":
    case "CXglyc":
    case "Condensation_DCL":
    case "Condensation_LCL":
    case "Condensation_Starter":
    case "Condensation_Dual":
    case "Heterocyclization":
        return "rgb(129,129,247)"
    case "Epimerization":
        return "rgb(129,129,247)"
    case "NRPS-COM_Nterm":
    case "NRPS-COM_Cterm":
    case "PKS_Docking_Nterm":
    case "PKS_Docking_Cterm":
    case "Trans-AT_docking":
        return "rgb(128,128,245)"
    case "Thioesterase":
    case "TD":
        return "rgb(245,196,242)"
    case "PKS_KS":
        return "rgb(129,247,129)"
    case "PKS_AT":
        return "rgb(247,129,129)"
    case "PKS_KR":
        return "rgb(128,246,128)"
    case "PKS_DH":
    case "PKS_DH2":
    case "PKS_DHt":
        return "rgb(247,190,129)"
    case "PKS_ER":
        return "rgb(129,247,243)"
    case "Aminotran_1_2":
    case "Aminotran_3":
    case "Aminotran_4":
    case "Aminotran_5":
    case "Polyketide_cyc2":
    default:
        return "rgb(218,218,218)"
    }
}


domain_label_type = function(type) {
    switch (type) {
    case "AMP-binding":
    case "AOX":
        return "A"
    case "PCP":
    case "ACP":
    case "NRPS-COM_Nterm":
    case "NRPS-COM_Cterm":
    case "PKS_Docking_Nterm":
    case "PKS_Docking_Cterm":
    case "Trans-AT_docking":
    case "Aminotran_1_2":
    case "Aminotran_3":
    case "Aminotran_4":
    case "Aminotran_5":
    case "Polyketide_cyc2":
        return ""
    case "Cglyc":
    case "CXglyc":
    case "Condensation_DCL":
    case "Condensation_LCL":
    case "Condensation_Starter":
    case "Condensation_Dual":
    case "Heterocyclization":
        return "C"
    case "Epimerization":
        return "E"
    case "Thioesterase":
        return "TE"
    case "PKS_KS":
        return "KS"
    case "PKS_AT":
        return "AT"
    case "PKS_KR":
        return "KR"
    case "PKS_DH":
    case "PKS_DH2":
        return "DH"
    case "PKS_DHt":
        return "DHt"
    case "PKS_ER":
        return "ER"
    default:
        return type.split('_')[0];
    }
}


function chunk_string(str, length) {
    if (typeof(str) !== 'undefined') {
        var chuncks = str.match(new RegExp('.{1,' + length + '}', 'g'));
        return chuncks.join('\n');
    }
    return '';
}


/* ***************************************************************************
 *  All the following functions are intented to change the positions of ORFs *
 *  belonging to the minus strand.                                           *
 *  Principle example:                                                       *
 *     0+ 1- 2- 3- 4+ -> 0+ 3- 2- 1- 4+                                      *
 *****************************************************************************/

function get_orf_strand(geneclusters, cluster_name, orf_name) {
    var cluster_name = cluster_name.replace(/-details/i, '');
    var orfs = geneclusters[cluster_name]['orfs'];
    for (var orf_idx in orfs) {
        if (orfs[orf_idx].locus_tag == orf_name) {
            return orfs[orf_idx].strand
        }
    }
}


function get_orfs_idx(cluster_name, details_data) {
    /*
     *  For a given cluster_name, find the index of
     *  ORFS in details_data structure
     */
    var orfs_idx = {}
    var orfs = details_data[cluster_name]['orfs'];
    for (var idx = 0; idx < orfs.length; idx++) {
        var orf_id = orfs[idx]['id'];
        orfs_idx[orf_id] = idx;
    }
    return orfs_idx;
}


function remove_values_if_notin(array1, array2) {
    /*
     *  Remove values from array1 if not in array2 -- keep order
     */
    var new_array = [];
    for (var idx = 0; idx < array1.length; idx++) {
        var elt = array1[idx];
        var index = array2.indexOf(elt);
        if ( index !== -1) {
            new_array.push(array1[idx]);
        }
    }
    return new_array;
}


function find_consecutive_minus_strand_orfs(cluster_name, geneclusters, details_data) {
    /* For a given cluster_name, find all consecutives
     *  ORFs belonging to minus strand
     *  Example:
     *  0+ 1- 2+ 3- 4- 5+ 6- 7- 8- --> [[3-, 4-], [6-, 7-, 8-]]
    */

    var cluster = geneclusters[cluster_name]
    //get orfs & index from details_data structure
    var orfs2idx = get_orfs_idx(cluster_name, details_data);
    var orfs_names = Object.keys(orfs2idx);

    var consecutive_minus_strand_orfs = [];
    var orfs = cluster['orfs'];
    var idx = 0;
    while (idx < orfs.length) {
        var minus_orfs = [];

        while ((idx < orfs.length) && (orfs[idx]['strand'] < 0)) {
            var orf_name = orfs[idx]['locus_tag'];
            minus_orfs.push(orf_name);
            idx++;
        }
        // remove ORFs no present in details_data (with no domains)
        minus_orfs = remove_values_if_notin(minus_orfs, orfs_names);

        // isolated ORFs belonging to minus strand are ignored
        if (minus_orfs.length > 1) {
            consecutive_minus_strand_orfs.push(minus_orfs);
        }
        idx++;
    }
    return consecutive_minus_strand_orfs;
}


function reverse_array(arr, i, j) {
    /*  Reverse inplace the chunk [i,j] of the
     *  array
     */
    if (i > j) { // ensure i <= j
        var tmp = i;
        i = j;
        j = tmp;
    }
    while (i < j) {
        var tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
        i++;
        j--;
    }
}


function is_sorted(arr) {
    for(var i = 0; i < (arr.length - 1); ++i) {
        if(arr[i] > arr[i+1]) {
            return false;
        }
    }
    return true;
}


function reverse_consecutive_minus_orfs(cluster_name, geneclusters, details_data) {
    /* For a given cluster_name, actually reverse consecutives
     *  ORFs belonging to minus strand
     *  Example:
     *  genecluster data: 0+ 1- 2+ 3- 4- 5- 6+ 7- 8- 9-
     *  get_orfs_idx: { 0: 0,
                        1: 1,
                        2: 2,
                        3: 3,
                        5: 4,
                        6: 5,
                        7: 6,
                        8: 7,
                        9: 8}
        (orf 4 in gencluster is no present in details_data)
     * consecutive_orfs: [[3-, 5-], [7-, 8-, 9-]]
     * details data before swap: 0+ 1- 2+ 3- 5- 6+ 7- 8- 9-
     * details data after swap: 0+ 1- 2+ 5- 3- 6+ 9- 8- 7-
     */

    var consecutive_orfs = find_consecutive_minus_strand_orfs(cluster_name, geneclusters, details_data);
    var orfs2idx = get_orfs_idx(cluster_name, details_data);
    var cluster = details_data[cluster_name];

    for (var idx = 0; idx < consecutive_orfs.length; idx++) {
        var chunk_names = consecutive_orfs[idx];
        var chunk_idx = [];
        for (var idx2 = 0; idx2 < chunk_names.length; idx2++) {
            var current_orf_name = chunk_names[idx2];
            var current_orf_idx = orfs2idx[current_orf_name];
            chunk_idx.push(current_orf_idx);
        }
        var start = Math.min(...chunk_idx);
        var end = Math.max(...chunk_idx);
        reverse_array(cluster['orfs'], start, end)

    }
}

function check_prerequisites(cluster_name, geneclusters, details_data) {
    /*
      To reverse minus orfs, some prerequistes are needed.
      First, geneclusters orfs must be order by start position
      Second, details_data orfs must follow the same order.
    */

    var start_positions = [];
    var geneclusters_orfs_names = []
    var orfs = geneclusters[cluster_name]['orfs'];
    for (var idx = 0; idx < orfs.length; idx++) {
        start_positions.push(orfs[idx].start);
        geneclusters_orfs_names.push(orfs[idx].locus_tag);
    }
    // assert orfs are ordered by start positions
    if (! is_sorted(start_positions)) {
        return false;
    }

    var orfs2idx = get_orfs_idx(cluster_name, details_data);
    var details_data_orfs_names = Object.keys(orfs2idx);
    geneclusters_orfs_names = remove_values_if_notin(geneclusters_orfs_names, details_data_orfs_names)
    if (geneclusters_orfs_names.toString() != details_data_orfs_names.toString())
    {
        return false;
    }
    return true;
}


$(document).ready(function() {
    $('#btn-direction').click(function() {

        //check prerequisites before
        for (var cluster_name in details_data) {
            var ok = check_prerequisites(cluster_name, geneclusters, details_data);
            // obviously, in swap mode, the orfs order for details_data is different from
            // order of geneclusters
            if ((! window.swapped) && (!ok)) {
                alert('This data did not meet the prerequisites to allow the swapping');
                return
            }
        }

        // Apply the reverse transformation on details_data
        // for all clusters & update the view
        for (var cluster_name in details_data) {
            reverse_consecutive_minus_orfs(cluster_name, geneclusters, details_data);
        }

        // if a cluster is already drawn, redraw it.
        if ( window.active_cluster != null ) {
            cluster_draw_svgs(details_data[window.active_cluster], 100, 700);
        }

	// if a domain is already drawn, deselect it.
	if ( window.active_cluster != null ) {
	    var lower = d3.select('#lower');
	    lower.selectAll('div').remove();
	    window.active_domain = null;
	}

        $( this ).toggleClass( "btn-danger" );
        $( this ).toggleClass( "btn-success" );

        if (window.swapped) {
            window.swapped = false;
        }
        else {
            window.swapped = true;
        }
    });
});



/* *************************************************************************
 *  All the following functions are intented build a fasta header foreach  *
 *  domains                                                                *
 ***************************************************************************/

class DefaultDict {
    /*
     * javascript equivalent of python defaultdict
     * https://stackoverflow.com/a/44622467
     */
  constructor(defaultInit) {
    return new Proxy({}, {
      get: (target, name) => name in target ?
        target[name] :
        (target[name] = typeof defaultInit === 'function' ?
          new defaultInit().valueOf() :
          defaultInit)
    })
  }
}


function sanitaze(field) {
    return field.replace(/[\W_]+/g, "-");
}

function complete_with_header(details_data) {
    /*
     * Build the fasta header for all domains
     * A/C/E domains will have a supplementary info:
     * An alternative domain denomination.
     * C A C A -> C1 A1 C2 A2 (count domains by type)
     */
    for (var cluster_name in details_data) {
	var orfs = details_data[cluster_name]['orfs'];
	for (var i = 0; i < orfs.length; i++) {
	    const counts = new DefaultDict(Number)
	    var domains = orfs[i]['domains'];
	    for (var j = 0; j < domains.length; j++) {
		var label = domain_label_type(domains[j]['type']);
		var mandatory_header = ['>' + cluster_name,
					sanitaze(orfs[i].id),
					j + 1]
		if (label == 'A' || label == 'C' || label == 'E') {
		    counts[label]++;
		    mandatory_header.push(label + counts[label]);

		}

		var optional_header = ['run_id=' + sanitaze(run_id),
				       'seq_name=' + sanitaze(retrieve_sequence_name(cluster_name, geneclusters)),
				       'type=' + sanitaze(domains[j]['type'])]

		domains[j].header = mandatory_header.join('_') + ' ' + optional_header.join(';');
	    }
	}
    }
}
