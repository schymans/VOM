function vom_init_vegpar(node, onCompleteCallback)
{
   var myobj = { label: "waterbalance_init", id: "waterbalance_init", href: "waterbalance_init.html", target:"basefrm" };
   var tmpNode24 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
