function watbal__f90(node, onCompleteCallback)
{
   var myobj = { label: "waterbalance", id: "waterbalance", href: "waterbalance.html", target:"basefrm" };
   var tmpNode66 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode66.isLeaf = true; 
   var myobj = { label: "waterbalance_init", id: "waterbalance_init", href: "waterbalance_init.html", target:"basefrm" };
   var tmpNode67 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode67.isLeaf = true; 
   var myobj = { label: "waterbalance_fluxes", id: "waterbalance_fluxes", href: "waterbalance_fluxes.html", target:"basefrm" };
   var tmpNode68 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode68.isLeaf = true; 
   var myobj = { label: "waterbalance_timestep", id: "waterbalance_timestep", href: "waterbalance_timestep.html", target:"basefrm" };
   var tmpNode69 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode69.isLeaf = true; 
   var myobj = { label: "waterbalance_update_state", id: "waterbalance_update_state", href: "waterbalance_update_state.html", target:"basefrm" };
   var tmpNode70 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode70.isLeaf = true; 
   var myobj = { label: "waterbalance_diag", id: "waterbalance_diag", href: "waterbalance_diag.html", target:"basefrm" };
   var tmpNode71 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode71.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
