function random_samples(node, onCompleteCallback)
{
   var myobj = { label: "read_shufflepar", id: "read_shufflepar", href: "read_shufflepar.html", target:"basefrm" };
   var tmpNode166 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode166.isLeaf = true; 
   var myobj = { label: "read_shufflevar", id: "read_shufflevar", href: "read_shufflevar.html", target:"basefrm" };
   var tmpNode167 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "transpmodel_init_once", id: "transpmodel_init_once", href: "transpmodel_init_once.html", target:"basefrm" };
   var tmpNode168 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "transpmodel", id: "transpmodel", href: "transpmodel.html", target:"basefrm" };
   var tmpNode179 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
