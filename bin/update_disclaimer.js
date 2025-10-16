#!/usr/bin/env node

const {EOL} = require('os');
const fs = require('fs');
const path = require('path');
const showdown  = require('showdown');

const filename = path.basename(__filename);
const docDir = path.resolve(filename, '../docs');
const hbsDir = path.resolve(filename, '../src/main/resources/org/pharmgkb/pharmcat/reporter')


try {
  let data = fs.readFileSync(path.resolve(docDir, 'Disclaimers.md'), 'utf8');
  data = data.replace(/---.*?---/s, '');

  const converter = new showdown.Converter();
  // don't use GitHub flavor because that adds incorrect line breaks (<br />)
  converter.setOption('ghCompatibleHeaderId', true);
  converter.setOption('headerLevelStart', 2);
  converter.setOption('simplifiedAutoLink', true);
  converter.setOption('ghCodeBlocks', true);
  converter.setOption('tables', true);
  converter.setOption('disableForced4SpacesIndentedSublists', true);

  const html = '<section id="disclaimer">' + EOL +
      converter.makeHtml(data)
          .replace('>Disclaimers and Other Information</h2>', '>Section IV: Disclaimers and Other Information</h2>')
      + EOL +
      '</section>' + EOL;
  fs.writeFileSync(path.resolve(hbsDir, 'disclaimers.hbs'), html);

} catch (err) {
  console.error(err);
}
