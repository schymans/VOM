function transpmodel_init_once(node, onCompleteCallback)
{
   var myobj = { label: "vom_read_input", id: "vom_read_input", href: "vom_read_input.html", target:"basefrm" };
   var tmpNode6 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_alloc", id: "vom_alloc", href: "vom_alloc.html", target:"basefrm" };
   var tmpNode9 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_open_output", id: "vom_open_output", href: "vom_open_output.html", target:"basefrm" };
   var tmpNode10 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode10.isLeaf = true; 
   var myobj = { label: "vom_get_soilprofile", id: "vom_get_soilprofile", href: "vom_get_soilprofile.html", target:"basefrm" };
   var tmpNode11 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_get_hourly_clim", id: "vom_get_hourly_clim", href: "vom_get_hourly_clim.html", target:"basefrm" };
   var tmpNode12 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "vom_get_perc_cov", id: "vom_get_perc_cov", href: "vom_get_perc_cov.html", target:"basefrm" };
   var tmpNode14 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
