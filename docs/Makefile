setup:
	bundle install --path vendor/bundle

serve:
	bundle exec jekyll serve --incremental

build:
	bundle exec jekyll build

clean:
	rm -rf _site

prototype: clean build
	surge _site
