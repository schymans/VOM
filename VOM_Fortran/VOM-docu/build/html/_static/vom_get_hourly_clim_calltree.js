function vom_get_hourly_clim(node, onCompleteCallback)
{
   var myobj = { label: "vom_calc_derived", id: "vom_calc_derived", href: "vom_calc_derived.html", target:"basefrm" };
   var tmpNode9 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
