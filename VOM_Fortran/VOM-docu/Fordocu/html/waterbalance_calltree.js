function waterbalance(node, onCompleteCallback)
{
   var myobj = { label: "waterbalance_fluxes", id: "waterbalance_fluxes", href: "waterbalance_fluxes.html", target:"basefrm" };
   var tmpNode35 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "waterbalance_timestep", id: "waterbalance_timestep", href: "waterbalance_timestep.html", target:"basefrm" };
   var tmpNode36 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "waterbalance_update_state", id: "waterbalance_update_state", href: "waterbalance_update_state.html", target:"basefrm" };
   var tmpNode37 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "waterbalance_diag", id: "waterbalance_diag", href: "waterbalance_diag.html", target:"basefrm" };
   var tmpNode38 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode38.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
